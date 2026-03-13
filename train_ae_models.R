# train_ae_models.R
# Builds autoencoder-based predictors per family & property, saving:
#   - <Property>.h5         (Keras model)
#   - <Property>_scaler.rds (list(mean, sd, method))
#   - <Property>_scaler.json
#
# Families:
#   ASD_DRY, NeoSpectra_DRY, NaturaSpec_DRY,
#   Agnostic_DRY (ASD + DRY subsets),
#   Agnostic_Moisture (ASD + all Neo/Natura moisture)
#
# Inputs (as provided by you):
#   data/VNIR to model/Scans_Track_NeoSpectra.xlsx (sheet="Scans_NeoSpectra")
#   data/VNIR to model/Scans_Track_NaturaSpec.xlsx (sheet="Scans_NaturaSpec")
#   data/VNIR to model/Scans_Track_ASD.xlsx       (sheet="Scans_ASD")
#   data/cornell_soil/Merged.xlsx
#
# Notes:
# - We align spectra to the target wavegrid per family (linear interpolation).
# - Optional SNV/SG preprocessing can be enabled to match app behavior.
# - Target y is z-scaled; scaler is saved to back-transform in Shiny.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(jsonlite)
  library(data.table)
  library(keras)
})

# -------------------- Configuration -------------------------------------------

soil_properties <- c(
  "soil_texture_sand","soil_texture_silt","soil_texture_clay",
  "organic_matter","soc","total_c","total_n",
  "active_carbon","ph","p","k",
  "mg","fe","mn","zn",
  "al","Ca","Cu","S",
  "B","pred_soil_protein","respiration","bd_ws"
)

# Wavegrids per sensor/family
WG_ASD_NATURA <- 350:2500
WG_NEO        <- 1350:2500
WG_INTERSECT  <- 1350:2500

split_idx <- function(n, seed = 42, p_train = 0.6, p_cal = 0.2, p_test = 0.2) {
  stopifnot(abs(p_train + p_cal + p_test - 1) < 1e-9)
  set.seed(seed); idx <- sample.int(n)
  n_tr <- floor(p_train * n); n_ca <- floor(p_cal * n)
  tr <- idx[seq_len(n_tr)]; ca <- idx[seq_len(n_ca) + n_tr]; te <- idx[-c(tr, ca)]
  list(train = tr, calib = ca, test = te)
}

metrics_from_y <- function(y_true, y_pred) {
  rmse <- sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
  sst  <- sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
  sse  <- sum((y_true - y_pred)^2, na.rm = TRUE)
  r2   <- if (sst > 0) 1 - sse/sst else NA_real_
  rpiq <- (IQR(y_true, na.rm = TRUE)) / rmse
  list(RMSE = rmse, R2 = r2, RPIQ = rpiq)
}

feature_minmax <- function(M) {
  data.frame(
    band = colnames(M),
    min = apply(M, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))),
    max = apply(M, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))),
    check.names = FALSE
  )
}

save_metrics_json <- function(path_base, metrics, conf_q, feat_range, latent) {
  obj <- list(
    metrics = metrics,                                  # RMSE, R2, RPIQ (test)
    conformal = list(q90 = conf_q[["q90"]], q95 = conf_q[["q95"]]),  # abs residual quantiles
    feature_range = feat_range,                         # post-preprocess min/max per band
    latent = latent                                     # μ, Σ, df, thr95
  )
  jsonlite::write_json(obj, paste0(path_base, "_metrics.json"), auto_unbox = TRUE, pretty = TRUE)
}

feature_minmax <- function(M) {
  data.frame(
    band = colnames(M),
    min = apply(M, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))),
    max = apply(M, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))),
    check.names = FALSE
  )
}

# Preprocessing to apply to X before training (must match app if enabled there)
# Absorbance, then SG with first derivative
PREPROCESS_STEPS <- c("ABSORBANCE","SG(11,2,1)")

EPOCHS      <- 30        # start modest; bump after confirming pipeline
BATCH_SIZE  <- 32
VAL_SPLIT   <- 0.2
PATIENCE_ES <- 8 # 8
PATIENCE_RL <- 4 # 4

# -------------------- Model builder (your AE) ---------------------------------

build_ae_asym <- function(d_in) {
  inp <- layer_input(shape = d_in)
  lat <- inp %>% 
    # layer_dense(256, activation="relu") %>% 
    layer_dense(128, activation="relu") %>% 
    layer_dense(64,  activation="relu") %>% 
    layer_dense(16,  activation="relu", name="latent")
  rec  <- lat %>% 
    layer_dense(32, activation="relu") %>% 
    layer_dense(d_in, activation="linear", name="reconstruction")
  pred <- lat %>% 
    layer_dense(64, activation="relu") %>% 
    layer_dropout(0.05) %>% 
    layer_dense(1, activation="linear", name="prediction")
  model <- keras_model(inputs = inp, outputs = list(rec, pred))
  model %>% compile(
    optimizer    = "adam",
    loss         = list("mse","mse"),
    loss_weights = list(0.3,0.3)
  )
}

# -------------------- Spectral utils (align + preprocess) ---------------------

get_wavelengths <- function(df, id_col = "Soil_ID") {
  cols <- setdiff(names(df), id_col)
  cols <- as.character(cols)
  wl   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", cols)))
  idx  <- !is.na(wl)
  list(wl = wl[idx], cols = cols[idx])
}

resample_to_grid <- function(M, src_wl, target_wl) {
  # M: numeric matrix (rows = samples, cols = src_wl order)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = length(target_wl))
  for (i in seq_len(nrow(M))) {
    out[i, ] <- approx(x = src_wl, y = M[i, ], xout = target_wl, rule = 1, ties = mean)$y
  }
  colnames(out) <- as.character(target_wl)
  out
}

# ----- Reflectance -> Absorbance ---------------------------------
# Converts R in [0,1] (or 0-100%) to absorbance A = -log(R)  (i.e., ln(1/R))
# Set base10 = TRUE if you prefer A = -log10(R).
reflectance_to_absorbance <- function(M, base10 = FALSE) {
  Mnum <- as.matrix(M)
  # If values look like percentages, scale to [0,1]
  if (is.finite(max(Mnum, na.rm = TRUE)) && max(Mnum, na.rm = TRUE) > 2) {
    Mnum <- Mnum / 100
  }
  # Avoid log(0)
  eps <- 1e-6
  Mnum[Mnum < eps] <- eps
  if (base10) {
    # A = log10(1/R) = -log10(R)
    return(-log10(Mnum))
  } else {
    # A = ln(1/R) = -ln(R)
    return(-log(Mnum))
  }
}

# ----- SG parsing now supports derivative d (e.g., "SG(11,2,1)") ---
parse_sg <- function(st) {
  inside <- substr(st, 4, nchar(st) - 1)   # "m,p" or "m,p,d"
  parts  <- strsplit(inside, ",", fixed = TRUE)[[1]]
  parts  <- trimws(parts)
  m <- suppressWarnings(as.integer(parts[1])); if (is.na(m)) m <- 11
  p <- suppressWarnings(as.integer(parts[2])); if (is.na(p)) p <- 2
  d <- if (length(parts) >= 3) suppressWarnings(as.integer(parts[3])) else 0
  if (is.na(d)) d <- 0
  list(m = m, p = p, d = d)
}

apply_sg_matrix <- function(M, m, p, d = 0) {
  ncols <- ncol(M)
  pars  <- sg_safe_params(m, p, ncols)
  if (is.null(pars)) return(M)
  m <- pars$m; p <- pars$p
  res <- t(apply(M, 1, function(r) {
    if (all(is.na(r))) return(r)
    # prospectr::savitzkyGolay: w = derivative order
    tryCatch(prospectr::savitzkyGolay(r, m = m, p = p, w = d),
             error = function(e) r)
  }))
  colnames(res) <- colnames(M)
  res
}

# ----- Pipeline: supports ABSORBANCE and SG(…, …, d) --------------
apply_pipeline <- function(M, steps, absorbance_base10 = FALSE) {
  out <- M
  if (length(steps) == 0) return(out)
  for (st in steps) {
    if (identical(st, "ABSORBANCE")) {
      out <- reflectance_to_absorbance(out, base10 = absorbance_base10)
    } else if (startsWith(st, "SG(")) {
      sg <- parse_sg(st)
      out <- apply_sg_matrix(out, m = sg$m, p = sg$p, d = sg$d)
    } else if (identical(st, "SNV")) {
      out <- apply_snv(out)
    }
  }
  out
}

sg_safe_params <- function(m, p, ncols) {
  if (!is.finite(m) || m < 3) m <- 3
  if (!is.finite(p) || p < 0) p <- 2
  if (m %% 2 == 0) m <- m + 1
  min_req <- p + 1
  if (m < min_req) m <- min_req + ifelse(min_req %% 2 == 0, 1, 0)
  if (ncols < 3) return(NULL)
  if (m > ncols) m <- ncols - ifelse(ncols %% 2 == 0, 1, 0)
  if (m < 3) return(NULL)
  list(m = m, p = p)
}


# -------------------- IO helpers ----------------------------------------------

dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

save_scaler <- function(mean_y, sd_y, path_base) {
  scaler <- list(mean = as.numeric(mean_y), sd = as.numeric(sd_y), method = "z")
  # RDS
  saveRDS(scaler, file = paste0(path_base, "_scaler.rds"))
  # JSON
  write_json(scaler, path = paste0(path_base, "_scaler.json"), auto_unbox = TRUE, pretty = TRUE)
}

# -------------------- Join & build X/Y ----------------------------------------

# Expects Soil_ID in all scan frames and soil table
# Will align spectra to target wavegrid, apply preprocessing, and return (X, y)
prepare_xy <- function(scans_df, soil_df, target_wg, prop, preprocess_steps = PREPROCESS_STEPS) {
  
  df <- scans_df %>%
    inner_join(soil_df, by = "Soil_ID")
  
  # extract spectra
  wlinfo <- get_wavelengths(df, id_col = "Soil_ID")
  stopifnot(length(wlinfo$wl) > 0)
  X_src <- as.matrix(df[, wlinfo$cols, drop = FALSE])
  
  # resample to target grid
  X <- resample_to_grid(X_src, src_wl = wlinfo$wl, target_wl = target_wg)
  
  # preprocess (optional)
  if (length(preprocess_steps) > 0) {
    X <- apply_pipeline(X, preprocess_steps)
  }
  
  # target
  y <- df[[prop]]
  keep <- is.finite(y)
  list(
    X = X[keep, , drop = FALSE],
    y = as.numeric(y[keep]),
    ids = df$Soil_ID[keep]
  )
}

# -------------------- Train one property model --------------------------------

train_property_model <- function(X, y, family_id, prop,
                                 epochs = EPOCHS, batch_size = BATCH_SIZE,
                                 val_split = VAL_SPLIT) {
  n <- nrow(X)
  if (length(y) < 30 || n < 30) {
    message("Skip ", family_id, " / ", prop, ": need ~30+, have ", n)
    return(invisible(NULL))
  }
  
  # scale target
  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- sd(y, na.rm = TRUE); if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z    <- as.numeric((y - mean_y) / sd_y)
  
  # split
  sp <- split_idx(n, seed = 42)
  Xtr <- X[sp$train, , drop = FALSE];  ytr_z <- y_z[sp$train]
  Xca <- X[sp$calib, , drop = FALSE];  yca   <- y[sp$calib]   # original units
  Xte <- X[sp$test,  , drop = FALSE];  yte   <- y[sp$test]    # original units
  
  # build & train
  d_in <- ncol(Xtr)
  mdl <- build_ae_asym(d_in)
  cb <- list(
    callback_early_stopping(monitor = "val_loss", patience = PATIENCE_ES, restore_best_weights = TRUE),
    callback_reduce_lr_on_plateau(monitor = "val_loss", patience = PATIENCE_RL, factor = 0.5, min_lr = 1e-5)
  )
  mdl %>% fit(
    x = as.matrix(Xtr),
    y = list(as.matrix(Xtr), ytr_z),
    epochs = epochs, batch_size = batch_size, validation_split = val_split,
    callbacks = cb, verbose = 0
  )
  
  # save model + scaler
  outdir <- file.path("models", family_id, "models")
  dir_create(outdir)
  base <- file.path(outdir, prop)
  keras::save_model_hdf5(mdl, filepath = paste0(base, ".h5"))
  save_scaler(mean_y, sd_y, path_base = base)
  
  # ---------- Conformal q90/q95 on calibration (orig units)
  pred_ca <- predict(mdl, as.matrix(Xca), verbose = 0)
  yhat_ca_z <- if (is.list(pred_ca)) {
    if (!is.null(names(pred_ca)) && "prediction" %in% names(pred_ca)) pred_ca[["prediction"]] else pred_ca[[length(pred_ca)]]
  } else pred_ca
  yhat_ca <- as.numeric(yhat_ca_z) * sd_y + mean_y
  res_ca  <- abs(yca - yhat_ca)
  q90 <- as.numeric(quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))
  
  # ---------- Test metrics (orig units)
  pred_te <- predict(mdl, as.matrix(Xte), verbose = 0)
  yhat_te_z <- if (is.list(pred_te)) {
    if (!is.null(names(pred_te)) && "prediction" %in% names(pred_te)) pred_te[["prediction"]] else pred_te[[length(pred_te)]]
  } else pred_te
  yhat_te <- as.numeric(yhat_te_z) * sd_y + mean_y
  mets <- metrics_from_y(y_true = yte, y_pred = yhat_te)
  
  # ---------- Latent stats on TRAIN (for applicability)
  encoder <- keras::keras_model(inputs = mdl$input,
                                outputs = mdl$get_layer("latent")$output)
  Ztr <- predict(encoder, as.matrix(Xtr), verbose = 0)  # n_tr x d
  mu_z <- colMeans(Ztr, na.rm = TRUE)
  Sigma_z <- stats::cov(Ztr, use = "pairwise.complete.obs")
  # small ridge to ensure PD
  eps <- 1e-6
  Sigma_z <- as.matrix(Sigma_z + diag(eps, ncol(Sigma_z)))
  d_lat <- length(mu_z)
  thr95 <- stats::qchisq(0.95, df = d_lat)
  
  latent <- list(mu = as.numeric(mu_z),
                 Sigma = unname(asplit(Sigma_z, 1L)),  # rows as list for compact JSON
                 df = d_lat, thr95 = thr95)
  
  # ---------- Feature min/max on full X (post-preprocess)
  fr <- feature_minmax(X)
  
  # save metrics JSON
  save_metrics_json(path_base = base,
                    metrics = mets,
                    conf_q  = list(q90 = q90, q95 = q95),
                    feat_range = fr,
                    latent = latent)
  
  invisible(TRUE)
}

# -------------------- Build a family (loop props) ------------------------------

build_family <- function(scans_df, soil_df, family_id, target_wg, soil_props) {
  message("=== Building family: ", family_id, " | grid ", min(target_wg), "-", max(target_wg), " (n=", length(target_wg), ") ===")
  for (prop in soil_props) {
    message(".. prep ", prop)
    xy <- prepare_xy(scans_df, soil_df, target_wg = target_wg, prop = prop)
    if (is.null(xy$X) || nrow(xy$X) == 0) {
      message("   (no usable rows)"); next
    }
    message(".. train ", prop, " (n=", nrow(xy$X), ", p=", ncol(xy$X), ")")
    try(train_property_model(xy$X, xy$y, family_id, prop), silent = FALSE)
  }
  message("=== Done: ", family_id, " ===")
}

# -------------------- Load data ------------------------------------------------

neospectra <- readxl::read_xlsx("Z:/Modeling/Climate_Smart_Modeling/data/VNIR to model/Scans_Track_NeoSpectra.xlsx", sheet = "Scans_NeoSpectra")
naturaspec <- readxl::read_xlsx("Z:/Modeling/Climate_Smart_Modeling/data/VNIR to model/Scans_Track_NaturaSpec.xlsx", sheet = "Scans_NaturaSpec")
asd        <- readxl::read_xlsx("Z:/Modeling/Climate_Smart_Modeling/data/VNIR to model/Scans_Track_ASD.xlsx", sheet = "Scans_ASD")
soil       <- readxl::read_xlsx("Z:/Modeling/Climate_Smart_Modeling/data/cornell_soil/Merged.xlsx")

# If moisture column names differ, adjust these predicates
is_dry <- function(df)   "DRY"  %in% names(df) && isTRUE(all(df$DRY  == 1)) # placeholder
is_1ML <- function(df)   "1ML"  %in% names(df) && isTRUE(all(df$`1ML` == 1))
is_3ML <- function(df)   "3ML"  %in% names(df) && isTRUE(all(df$`3ML` == 1))
# More robust: if you have a 'Moisture' factor/label column:
# filter(df, Moisture == "DRY") etc.

# If your moisture is encoded as a column "Moisture" with values "DRY","1ML","3ML":
neo_dry  <- if ("Lab_Treatment" %in% names(neospectra)) dplyr::filter(neospectra, Lab_Treatment == "DRY") else neospectra
nat_dry  <- if ("Lab_Treatment" %in% names(naturaspec)) dplyr::filter(naturaspec, Lab_Treatment == "DRY") else naturaspec

# -------------------- Build requested families --------------------------------

# ASD DRY
build_family(asd, soil, "ASD_DRY", WG_ASD_NATURA, soil_properties)

# NeoSpectra DRY
build_family(neo_dry, soil, "NeoSpectra_DRY", WG_NEO, soil_properties)

# NaturaSpec DRY
build_family(nat_dry, soil, "NaturaSpec_DRY", WG_ASD_NATURA, soil_properties)

# Agnostic DRY (ASD + DRY only from the others)
neo_dry$Field_Replicate <- as.numeric(neo_dry$Field_Replicate)
neo_dry$Lab_Replicate <- as.numeric(neo_dry$Lab_Replicate)

agnostic_dry <- dplyr::bind_rows(asd, neo_dry, nat_dry)
build_family(agnostic_dry, soil, "Agnostic_DRY", WG_INTERSECT, soil_properties)

# Agnostic Moisture (ASD + all Neo/Natura moisture levels pooled)
neospectra$Field_Replicate <- as.numeric(neospectra$Field_Replicate)
neospectra$Lab_Replicate <- as.numeric(neospectra$Lab_Replicate)

agnostic_all <- dplyr::bind_rows(asd, neospectra, naturaspec)

# Choose a wavegrid; if you want broadest coverage for ASD/NaturaSpec, use 350:2500
# (NeoSpectra will be interpolated below 1350 — OK for AE but mind extrapolation quality.)
build_family(agnostic_all, soil, "Agnostic_Moisture", WG_INTERSECT, soil_properties)

cat("\nAll families processed. Models + scalers saved under ./models/<FAMILY>/models/\n")
