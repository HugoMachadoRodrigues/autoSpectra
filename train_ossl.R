# train_ossl.R - Training orchestrator for autoSpectra (Safanelli 2025 faithful)
#
# Trains four OSSL model families:
#   OSSL_VisNIR        - masked soilVAE on 350-2500 nm
#   OSSL_MIR           - masked soilVAE on 600-4000 cm-1
#   OSSL_VisNIR_Cubist - Cubist + PCA(120) + SNV (Safanelli replication)
#   OSSL_MIR_Cubist    - Cubist + PCA(120) + SNV (Safanelli replication)
#
# Methodology per Safanelli et al. (2025) PLoS ONE 20(1):e0296545:
#   - SNV preprocessing (no derivatives)
#   - log1p target transform for skewed properties
#   - 10-fold external CV with refitting (outer loop)
#   - 5-fold inner CV for hyperparameter tuning (Cubist committees)
#   - Applicability domain via Q-statistic (Cubist) or Mahalanobis (VAE)
#   - Normalised conformal prediction at 68% confidence (Cubist)
#
# Outputs per (family, property):
#   models/<family_id>/models/<prop>.keras           (VAE) or
#   models/<family_id>/models/<prop>_cubist.rds      (Cubist)
#   models/<family_id>/models/<prop>_scaler.rds
#   models/<family_id>/models/<prop>_metrics.json
#   models/<family_id>/metrics_summary.json
#   models/<family_id>/cv10_summary.json  (10-fold CV aggregate metrics)
#
# Usage:
#   Rscript train_ossl.R                            # all 4 families, full
#   Rscript train_ossl.R OSSL_MIR                   # one family only
#   Rscript train_ossl.R --only-vae                 # both VAE families
#   Rscript train_ossl.R --only-cubist              # both Cubist families
#   Rscript train_ossl.R --cv10                     # 10-fold CV refit mode
#   Rscript train_ossl.R --quick                    # 3-fold, 20 epochs, 4 props
#   Rscript train_ossl.R --props oc,clay.tot,ph.h2o # selected properties

# Force CPU for VAE training - tensorflow-metal 1.2.0 + our masked_mse / Sampling
# layer combo produces NaN-collapsed predictions on Apple Silicon (verified
# 2026-04-26: VisNIR oc R^2 = -1.17 with metal vs +0.46 on CPU). Setting both
# Apple's metal-disable flag and the generic GPU-visibility env var.
Sys.setenv(DISABLE_TENSORFLOW_METAL = "1")
Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "2")

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) {
    if (requireNamespace("pkgload", quietly = TRUE)) {
      pkgload::load_all(".", quiet = TRUE)
    } else {
      for (f in list.files("R", pattern = "\\.R$", full.names = TRUE))
        source(f)
    }
  } else {
    library(autoSpectra)
  }
})

# ---- Argument parsing ----------------------------------------------------

args_raw <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) flag %in% args_raw
arg_val  <- function(flag, default = NULL) {
  i <- which(args_raw == flag)
  if (length(i) == 0 || i[1L] >= length(args_raw)) return(default)
  args_raw[i[1L] + 1L]
}

QUICK       <- has_flag("--quick")
CV10        <- has_flag("--cv10")
LODO        <- has_flag("--lodo")
ONLY_VAE    <- has_flag("--only-vae")
ONLY_CUBIST <- has_flag("--only-cubist")

PROPS_ARG <- arg_val("--props", NULL)
PROPS     <- if (!is.null(PROPS_ARG)) trimws(strsplit(PROPS_ARG, ",")[[1]]) else NULL

# --workers for parallelism (currently unused: keras + cubist are sequential
# per-property to avoid GPU contention; retained for future parallel support)
N_WORKERS <- suppressWarnings(as.integer(arg_val("--workers", "1")))
if (is.na(N_WORKERS)) N_WORKERS <- 1L

# Resolve family list from positional args (everything not starting with --
# and not the value of a known flag)
flag_values <- c(arg_val("--props"), arg_val("--workers"))
ALL_FAMILIES <- c("OSSL_VisNIR", "OSSL_MIR",
                  "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist")
POS <- args_raw[!grepl("^--", args_raw) & !(args_raw %in% flag_values)]
FAMILIES <- if (length(POS) > 0L) POS else ALL_FAMILIES
if (ONLY_VAE)    FAMILIES <- intersect(FAMILIES, c("OSSL_VisNIR", "OSSL_MIR"))
if (ONLY_CUBIST) FAMILIES <- intersect(FAMILIES,
                                       c("OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist"))

QUICK_PROPS <- c("oc", "clay.tot", "ph.h2o", "n.tot")

# ---- Configuration -------------------------------------------------------

K_FOLDS   <- if (QUICK) 3L  else 10L
EPOCHS    <- if (QUICK) 20L else 60L     # full: 60 (not 100) - ES typically triggers ~25-30
BATCH     <- if (QUICK) 64L else 128L    # bigger batch on GPU = faster per-epoch
PAT_ES    <- if (QUICK) 5L  else 12L     # full: 12 (not 20)
PAT_LR    <- if (QUICK) 3L  else 5L
LATENT    <- 32L          # tweaked: 16 -> 32 for more representational capacity
KL_BETA   <- 5e-4         # tweaked: 1e-3 -> 5e-4 (less KL pressure -> better recon/pred)
KL_WARMUP <- if (QUICK) 10L else 20L
N_AUG     <- if (QUICK) 2L  else 2L      # 2 aug copies (3 was too slow on full data)
MIN_N     <- 50L
OUT_DIR   <- "models"
CACHE_DIR <- autoSpectra::ossl_cache_dir()

message("=============================================================")
message("autoSpectra training - Safanelli 2025 faithful methodology")
message("=============================================================")
message("  Families      : ", paste(FAMILIES, collapse = ", "))
.props_lbl <- if (!is.null(PROPS)) paste(PROPS, collapse = ", ") else
              if (QUICK)          paste(QUICK_PROPS, collapse = ", ") else
                                  "all (34)"
.mode_lbl  <- if (CV10)  "10-fold CV refit" else
              if (QUICK) "QUICK"             else
                         "single train/calib/test split"
.fold_lbl  <- if (LODO)  "LODO (stratified by dataset)" else "random k-fold"
message("  Properties    : ", .props_lbl)
message("  Mode          : ", .mode_lbl)
message("  Fold strategy : ", .fold_lbl)

# ---- Per-family training -------------------------------------------------

.safe_nm <- function(family_id, mode = c("single", "cv10")) {
  mode <- match.arg(mode)
  dir_create(file.path(OUT_DIR, family_id))
  file.path(OUT_DIR, family_id,
            if (mode == "cv10") "cv10_summary.json" else "metrics_summary.json")
}

.collect_per_property_metrics <- function(family_id, properties) {
  rows <- list()
  for (prop in properties) {
    fp <- file.path(OUT_DIR, family_id, "models",
                    paste0(prop, "_metrics.json"))
    if (!file.exists(fp)) next
    info <- tryCatch(jsonlite::read_json(fp), error = function(e) NULL)
    if (is.null(info)) next
    m <- info$metrics
    rows[[prop]] <- list(
      property = prop,
      RMSE     = unlist(m$RMSE),
      Bias     = unlist(m$Bias),
      R2       = unlist(m$R2),
      RPIQ     = unlist(m$RPIQ),
      CCC      = unlist(m$CCC)
    )
  }
  rows
}

.summarise_family <- function(family_id, results, cv_summaries = NULL) {
  ok <- !vapply(results, is.null, logical(1))
  per_prop <- .collect_per_property_metrics(family_id, names(results)[ok])
  summary <- list(
    family_id     = family_id,
    n_properties  = length(results),
    n_successful  = sum(ok),
    n_failed      = sum(!ok),
    properties    = names(results),
    per_property  = per_prop
  )
  if (!is.null(cv_summaries) && length(cv_summaries) > 0L) {
    summary$cv10 <- cv_summaries
  }
  out_path <- .safe_nm(family_id, mode = if (CV10) "cv10" else "single")
  jsonlite::write_json(summary, out_path, auto_unbox = TRUE, pretty = TRUE)
  message("  Summary saved: ", out_path)
}

.prepare_preprocessed_matrix <- function(fam, joined_df) {
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)
  mask    <- mask_from_coverage(X_res)
  X_res[!is.finite(X_res)] <- 0
  X_proc  <- apply_pipeline(X_res, fam$preprocess, mask = mask)
  list(X = X_proc, mask = mask, wavegrid = fam$wavegrid)
}

.select_y <- function(joined_df, prop) {
  y <- suppressWarnings(as.numeric(joined_df[[prop]]))
  if (prop %in% ossl_l1_log_transform) {
    ok <- is.finite(y) & y > -1
    y[!ok] <- NA_real_
    y <- log1p(y)
  }
  y
}

.train_family_single <- function(family_id, props, joined_df) {
  fam  <- get_family(family_id)
  algo <- fam$algorithm %||% "vae"

  pp   <- .prepare_preprocessed_matrix(fam, joined_df)
  X    <- pp$X

  results <- list()
  for (prop in props) {
    y    <- .select_y(joined_df, prop)
    keep <- is.finite(y) & is.finite(rowSums(X))
    if (sum(keep) < MIN_N) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      next
    }
    message("  Training [", prop, "] (n=", sum(keep), ") ...")
    res <- try({
      if (identical(algo, "vae_masked")) {
        train_maskedVAE(X[keep, ], y[keep], family_id, prop,
                        wavegrid    = fam$wavegrid,
                        sensor_type = fam$sensor_type,
                        out_dir     = OUT_DIR,
                        epochs      = EPOCHS,
                        batch_size  = BATCH,
                        patience_es = PAT_ES,
                        patience_rl = PAT_LR,
                        latent_dim  = LATENT,
                        kl_beta     = KL_BETA,
                        kl_warmup_epochs = KL_WARMUP,
                        n_aug_copies = N_AUG,
                        min_n       = MIN_N)
      } else if (identical(algo, "cubist_pca")) {
        train_cubist_pca(X[keep, ], y[keep], family_id, prop,
                         out_dir = OUT_DIR, min_n = MIN_N)
      } else {
        train_soilVAE(X[keep, ], y[keep], family_id, prop,
                      out_dir = OUT_DIR,
                      epochs = EPOCHS, batch_size = BATCH,
                      patience_es = PAT_ES, patience_rl = PAT_LR,
                      latent_dim = LATENT, kl_beta = KL_BETA,
                      kl_warmup_epochs = KL_WARMUP, min_n = MIN_N)
      }
    }, silent = FALSE)
    results[[prop]] <- res
  }
  results
}

.train_family_cv10 <- function(family_id, props, joined_df) {
  fam  <- get_family(family_id)
  pp   <- .prepare_preprocessed_matrix(fam, joined_df)
  X    <- pp$X

  cv_summaries <- list()
  for (prop in props) {
    y    <- .select_y(joined_df, prop)
    keep <- is.finite(y) & is.finite(rowSums(X))
    if (sum(keep) < MIN_N) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      next
    }
    y_k  <- y[keep]; X_k <- X[keep, , drop = FALSE]

    folds <- if (LODO &&
                 "dataset.code_ascii_txt" %in% names(joined_df)) {
      ds <- joined_df[["dataset.code_ascii_txt"]][keep]
      stratified_kfold(ds, k = K_FOLDS, seed = 42L)
    } else {
      random_kfold(length(y_k), k = K_FOLDS, seed = 42L)
    }

    message("  10-fold CV [", prop, "] (n=", length(y_k), ") ...")
    cv_res <- try(
      cv_train_property(X_k, y_k, family_id, prop, folds, fam,
                        epochs = EPOCHS, batch_size = BATCH,
                        patience_es = PAT_ES, patience_rl = PAT_LR,
                        latent_dim = LATENT, kl_beta = KL_BETA,
                        kl_warmup_epochs = KL_WARMUP),
      silent = FALSE
    )
    cv_summaries[[prop]] <- if (inherits(cv_res, "try-error")) NULL else cv_res
    ms <- cv_summaries[[prop]]$mean_sd
    if (!is.null(ms)) {
      message(sprintf(
        "    RMSE=%.3f+-%.3f  R2=%.3f+-%.3f  RPIQ=%.2f+-%.2f  CCC=%.3f+-%.3f",
        ms$RMSE_mean, ms$RMSE_sd, ms$R2_mean, ms$R2_sd,
        ms$RPIQ_mean, ms$RPIQ_sd, ms$CCC_mean, ms$CCC_sd))
    }
  }
  cv_summaries
}

# ---- Main dispatch -------------------------------------------------------

# Resolve property list
target_props <- if (!is.null(PROPS)) PROPS else
                if (QUICK)           QUICK_PROPS else
                                     ossl_l1_properties

# Load OSSL data once per sensor_type (same data feeds VAE + Cubist variants)
joined_cache <- new.env(parent = emptyenv())
.load_joined <- function(sensor_type) {
  key <- sensor_type
  if (!exists(key, envir = joined_cache, inherits = FALSE)) {
    df <- ossl_prepare(sensor_type  = sensor_type,
                       cache_dir    = CACHE_DIR,
                       properties   = target_props,
                       download_if_missing = TRUE)
    assign(key, df, envir = joined_cache)
  }
  get(key, envir = joined_cache, inherits = FALSE)
}

for (fid in FAMILIES) {
  fam <- get_family(fid)
  message("")
  message("=============================================================")
  message("[", fid, "] ", fam$label)
  message("  Algorithm : ", fam$algorithm %||% "vae")
  message("  Pipeline  : ", paste(fam$preprocess, collapse = " -> "))
  message("=============================================================")

  joined_df <- .load_joined(fam$sensor_type)

  props <- intersect(target_props, fam$properties)

  if (CV10) {
    cv_summaries <- .train_family_cv10(fid, props, joined_df)
    # Mark "results" as done for whichever folds produced metrics
    results <- setNames(lapply(cv_summaries, function(s)
      if (is.null(s)) NULL else TRUE), names(cv_summaries))
    .summarise_family(fid, results, cv_summaries = cv_summaries)
  } else {
    results <- .train_family_single(fid, props, joined_df)
    .summarise_family(fid, results)
  }
}

message("")
message("=============================================================")
message("All requested families trained.")
message("Next step: Rscript upload_models.R  (package zips for Zenodo)")
message("=============================================================")
