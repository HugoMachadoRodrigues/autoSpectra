# train_ossl.R — Scientific training script for autoSpectra OSSL models
#
# Trains one soilVAE model per soil property × spectral domain using the full
# OSSL v1.2 corpus. Sensor-agnosticism is achieved through:
#   1. Multi-instrument training corpus (all OSSL contributing datasets)
#   2. Two-step SG first-derivative preprocessing (removes baseline +
#      multiplicative scatter — the two main sources of inter-instrument
#      spectral variation)
#   3. soilVAE information bottleneck (16D latent space forces encoder to
#      learn soil-chemistry features, discarding instrument artifacts)
#   4. Mahalanobis applicability domain in latent space
#
# Validation strategy:
#   10-fold cross-validation STRATIFIED BY OSSL CONTRIBUTING DATASET so that
#   each fold holds out a geographic/instrument cluster. This is strictly
#   harder than random CV and gives an honest estimate of cross-sensor
#   generalisation.
#
# Outputs (per family):
#   models/<family_id>/models/<prop>.h5           — soilVAE weights
#   models/<family_id>/models/<prop>_scaler.rds   — z-score parameters
#   models/<family_id>/models/<prop>_metrics.json — test metrics + latent stats
#   models/<family_id>/metrics_summary.json       — all properties, mean±SD
#
# Usage:
#   Rscript train_ossl.R                     # both families
#   Rscript train_ossl.R OSSL_VisNIR         # VisNIR only
#   Rscript train_ossl.R OSSL_MIR            # MIR only

suppressPackageStartupMessages({
  # Always load from local source so fixes take effect without re-installing.
  # Assumes the script is run from the package root (where DESCRIPTION lives).
  if (file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(autoSpectra)
  }
  library(keras3)
})

# ---- Configuration -------------------------------------------------------

args      <- commandArgs(trailingOnly = TRUE)
FAMILIES  <- if (length(args) > 0) args else c("OSSL_VisNIR", "OSSL_MIR")
K_FOLDS   <- 10L        # cross-validation folds (stratified by dataset)
EPOCHS    <- 80L        # max training epochs per fold
BATCH     <- 32L        # mini-batch size
PAT_ES    <- 10L        # early stopping patience
PAT_LR    <- 5L         # LR reduction patience
LATENT    <- 16L        # latent dimension
MIN_N     <- 50L        # minimum samples to attempt training
OUT_DIR   <- "models"
CACHE_DIR <- autoSpectra::ossl_cache_dir()

# ---- Helper: k-fold indices stratified by dataset ------------------------

#' @param dataset_col Character vector of dataset labels (one per sample)
#' @param k Number of folds
#' @param seed RNG seed
#' @return List of length k; each element is an integer vector of test indices
stratified_kfold <- function(dataset_col, k = 10L, seed = 42L) {
  set.seed(seed)
  folds <- vector("list", k)
  datasets <- unique(dataset_col)
  # Assign each dataset to a fold (round-robin to balance fold sizes)
  ds_order <- sample(datasets)
  ds_fold  <- setNames(((seq_along(ds_order) - 1L) %% k) + 1L, ds_order)
  for (fold in seq_len(k)) {
    in_fold  <- which(ds_fold[dataset_col] == fold)
    # If a dataset is too large to sit in a single fold, sub-sample within it
    folds[[fold]] <- in_fold
  }
  # Fallback: if a fold is empty (too few datasets), use random assignment
  empty <- which(vapply(folds, length, integer(1)) == 0L)
  if (length(empty) > 0) {
    message("  Note: ", length(empty), " empty fold(s) — ",
            "using random assignment instead of dataset stratification.")
    n   <- length(dataset_col)
    idx <- sample.int(n)
    cuts <- round(seq(0, n, length.out = k + 1))
    for (fold in empty)
      folds[[fold]] <- idx[(cuts[fold] + 1):cuts[fold + 1]]
  }
  folds
}

# ---- Helper: single-property CV training ---------------------------------

train_property_cv <- function(X, y, dataset_col, family_id, prop,
                               k = K_FOLDS, epochs = EPOCHS,
                               batch = BATCH, latent = LATENT,
                               pat_es = PAT_ES, pat_lr = PAT_LR,
                               out_dir = OUT_DIR) {
  n <- nrow(X)
  if (n < MIN_N) {
    message("    SKIP: need >= ", MIN_N, ", have ", n)
    return(NULL)
  }

  # Z-score scaler (fit on all data so final model is comparable)
  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  # K-fold CV — collect per-fold metrics
  folds      <- stratified_kfold(dataset_col, k = k)
  fold_rmse  <- numeric(k)
  fold_r2    <- numeric(k)
  fold_rpiq  <- numeric(k)

  for (fold in seq_len(k)) {
    te_idx  <- folds[[fold]]
    if (length(te_idx) == 0) {
      fold_rmse[fold] <- NA; fold_r2[fold] <- NA; fold_rpiq[fold] <- NA
      next
    }
    tr_idx <- setdiff(seq_len(n), te_idx)

    X_tr <- X[tr_idx, , drop = FALSE]; y_tr <- y_z[tr_idx]
    X_te <- X[te_idx, , drop = FALSE]; y_te <- y[te_idx]

    mdl_fold <- build_soilVAE(ncol(X_tr), latent_dim = latent)
    cbs <- list(
      keras3::callback_early_stopping(
        monitor = "val_loss", patience = pat_es,
        restore_best_weights = TRUE),
      keras3::callback_reduce_lr_on_plateau(
        monitor = "val_loss", patience = pat_lr,
        factor = 0.5, min_lr = 1e-5)
    )
    suppressMessages(
      mdl_fold |> keras3::fit(
        x = as.matrix(X_tr),
        y = list(as.matrix(X_tr), y_tr),
        epochs = epochs, batch_size = batch,
        validation_split = 0.15,
        callbacks = cbs, verbose = 0
      )
    )
    yhat_z      <- .extract_prediction(mdl_fold, X_te)
    yhat        <- yhat_z * sd_y + mean_y
    mets        <- metrics_from_y(y_te, yhat)
    fold_rmse[fold] <- mets$RMSE
    fold_r2[fold]   <- mets$R2
    fold_rpiq[fold] <- mets$RPIQ
    keras3::clear_session()   # free GPU memory between folds
    rm(mdl_fold); gc(verbose = FALSE)
  }

  cv <- list(
    RMSE_mean = mean(fold_rmse, na.rm = TRUE),
    RMSE_sd   = stats::sd(fold_rmse, na.rm = TRUE),
    R2_mean   = mean(fold_r2,   na.rm = TRUE),
    R2_sd     = stats::sd(fold_r2,   na.rm = TRUE),
    RPIQ_mean = mean(fold_rpiq, na.rm = TRUE),
    RPIQ_sd   = stats::sd(fold_rpiq, na.rm = TRUE),
    n_folds   = k,
    fold_RMSE = fold_rmse,
    fold_R2   = fold_r2,
    fold_RPIQ = fold_rpiq
  )

  # --- Final model: train on ALL data, split 80/10/10 for calib + latent ---
  sp     <- split_idx(n, seed = 42L, p_train = 0.8, p_cal = 0.1, p_test = 0.1)
  X_tr   <- X[sp$train, , drop = FALSE]; y_tr <- y_z[sp$train]
  X_ca   <- X[sp$calib, , drop = FALSE]; y_ca <- y[sp$calib]
  X_te   <- X[sp$test,  , drop = FALSE]; y_te <- y[sp$test]

  mdl_final <- build_soilVAE(ncol(X_tr), latent_dim = latent)
  cbs_final <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = pat_es,
      restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = pat_lr,
      factor = 0.5, min_lr = 1e-5)
  )
  suppressMessages(
    mdl_final |> keras3::fit(
      x = as.matrix(X_tr),
      y = list(as.matrix(X_tr), y_tr),
      epochs = epochs, batch_size = batch,
      validation_split = 0.15,
      callbacks = cbs_final, verbose = 0
    )
  )

  # Save model + scaler
  model_dir_prop <- file.path(out_dir, family_id, "models")
  dir_create(model_dir_prop)
  base <- file.path(model_dir_prop, prop)
  keras3::save_model(mdl_final, paste0(base, ".h5"))
  save_scaler(mean_y, sd_y, path_base = base)

  # Conformal calibration intervals (absolute residuals on calib set)
  yhat_ca <- .extract_prediction(mdl_final, X_ca) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(stats::quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(stats::quantile(res_ca, 0.95, na.rm = TRUE))

  # Test-set metrics (held-out 10%)
  yhat_te <- .extract_prediction(mdl_final, X_te) * sd_y + mean_y
  mets_te <- metrics_from_y(y_te, yhat_te)

  # Latent statistics for Mahalanobis applicability domain
  encoder <- keras3::keras_model(inputs  = mdl_final$input,
                                 outputs = mdl_final$get_layer("latent")$output)
  Z_tr  <- predict(encoder, as.matrix(X[sp$train, , drop = FALSE]), verbose = 0)
  mu_z  <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z <- stats::cov(Z_tr, use = "pairwise.complete.obs")
  Sig_z <- as.matrix(Sig_z + diag(1e-6, ncol(Sig_z)))
  d_lat <- length(mu_z)
  thr95 <- stats::qchisq(0.95, df = d_lat)

  latent_stats <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = d_lat,
    thr95 = thr95
  )

  save_metrics_json(base,
                    metrics    = mets_te,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = feature_minmax(X),
                    latent     = latent_stats)

  list(cv = cv, test = mets_te, n = n)
}

# ---- Main training loop --------------------------------------------------

for (fam_id in FAMILIES) {
  fam <- get_family(fam_id)
  message("\n", strrep("=", 60))
  message("Training: ", fam$label)
  message(strrep("=", 60))
  message("Grid     : ", length(fam$wavegrid), " bands (",
          min(fam$wavegrid), " - ", max(fam$wavegrid), ")")
  message("Pipeline : ", paste(fam$preprocess, collapse = " -> "))
  message("CV folds : ", K_FOLDS, " (stratified by OSSL dataset)")

  # ---- Load and prepare OSSL data ----------------------------------------
  message("\nLoading OSSL data ...")
  joined_df <- ossl_prepare(
    sensor_type         = fam$sensor_type,
    cache_dir           = CACHE_DIR,
    properties          = fam$properties,
    download_if_missing = TRUE
  )
  message("  Samples after join: ", nrow(joined_df))

  # Dataset column for stratification (use OSSL dataset_code if present)
  ds_col <- intersect(c("dataset.code_ascii_txt", "dataset_code", "location_id"),
                      names(joined_df))
  if (length(ds_col) == 0) {
    message("  WARNING: no dataset column found — using random CV folds")
    dataset_labels <- as.character(seq_len(nrow(joined_df)))
  } else {
    dataset_labels <- as.character(joined_df[[ds_col[1]]])
    n_ds <- length(unique(dataset_labels))
    message("  Datasets found: ", n_ds,
            " (", paste(head(unique(dataset_labels), 5), collapse = ", "),
            if (n_ds > 5) ", ..." else "", ")")
  }

  # Resample to canonical grid
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)

  # Preprocessing (removes inter-instrument baseline + scatter)
  message("  Preprocessing ...")
  X_proc <- apply_pipeline(X_res, fam$preprocess)

  # ---- Per-property training ---------------------------------------------
  summary_list <- list()

  for (prop in fam$properties) {
    y <- suppressWarnings(as.numeric(joined_df[[prop]]))
    keep <- is.finite(y) & is.finite(rowSums(X_proc))
    n_ok <- sum(keep)

    if (n_ok < MIN_N) {
      message(sprintf("  [%-20s]  SKIP  (n=%d)", prop, n_ok))
      next
    }

    message(sprintf("  [%-20s]  n=%-6d  training ...", prop, n_ok))

    result <- tryCatch(
      train_property_cv(
        X          = X_proc[keep, , drop = FALSE],
        y          = y[keep],
        dataset_col = dataset_labels[keep],
        family_id  = fam_id,
        prop       = prop
      ),
      error = function(e) {
        message("    ERROR: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(result)) {
      cv  <- result$cv
      message(sprintf(
        "    CV R2=%.3f\u00b1%.3f  RMSE=%.3f\u00b1%.3f  RPIQ=%.2f\u00b1%.2f",
        cv$R2_mean, cv$R2_sd,
        cv$RMSE_mean, cv$RMSE_sd,
        cv$RPIQ_mean, cv$RPIQ_sd
      ))
      summary_list[[prop]] <- list(
        property   = prop,
        label      = property_label(prop),
        n          = result$n,
        R2_mean    = cv$R2_mean,
        R2_sd      = cv$R2_sd,
        RMSE_mean  = cv$RMSE_mean,
        RMSE_sd    = cv$RMSE_sd,
        RPIQ_mean  = cv$RPIQ_mean,
        RPIQ_sd    = cv$RPIQ_sd
      )
    }
    keras3::clear_session()
    gc(verbose = FALSE)
  }

  # ---- Save family-level metrics summary ----------------------------------
  summary_path <- file.path(OUT_DIR, fam_id, "metrics_summary.json")
  dir_create(dirname(summary_path))
  jsonlite::write_json(list(
    family       = fam_id,
    ossl_version = fam$ossl_version,
    k_folds      = K_FOLDS,
    strategy     = "stratified by OSSL contributing dataset",
    properties   = summary_list
  ), summary_path, auto_unbox = TRUE, pretty = TRUE)

  message("\nSummary saved: ", summary_path)
  message("Done: ", fam_id, " (", length(summary_list), " properties trained)\n")
}

message("\n", strrep("=", 60))
message("All families complete. Run run_autoSpectra() to launch the app.")
message(strrep("=", 60))

autoSpectra::run_autoSpectra()
