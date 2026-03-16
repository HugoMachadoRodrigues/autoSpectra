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
#   Leave-One-Dataset-Out (LODO) cross-validation: each fold withholds all
#   samples from one OSSL contributing dataset (geographic/instrument cluster).
#   This matches the methodology of Safanelli et al. and is the most honest
#   estimate of cross-sensor generalisation in the soil spectroscopy literature.
#
# Metrics: RMSE, Bias, R², RPIQ (Bellon-Maurel 2010), CCC (Lin 1989)
#
# Outputs (per family):
#   models/<family_id>/models/<prop>.h5           — soilVAE weights
#   models/<family_id>/models/<prop>_scaler.rds   — z-score parameters
#   models/<family_id>/models/<prop>_metrics.json — test metrics + latent stats
#   models/<family_id>/metrics_summary.json       — all properties, mean±SD
#
# Usage:
#   Rscript train_ossl.R                      # both families, sequential
#   Rscript train_ossl.R OSSL_VisNIR          # VisNIR only
#   Rscript train_ossl.R --quick              # 3-fold, 20 epochs, 4 properties
#   Rscript train_ossl.R --parallel           # parallel property training
#   Rscript train_ossl.R --parallel --workers 4 OSSL_VisNIR

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(autoSpectra)
  }
  library(keras3)
})

# ---- Argument parsing -------------------------------------------------------

args_raw  <- commandArgs(trailingOnly = TRUE)

QUICK     <- "--quick"    %in% args_raw
PARALLEL  <- "--parallel" %in% args_raw

# --workers N  (default: physical cores, max 8)
.wi <- which(args_raw == "--workers")
N_WORKERS <- if (length(.wi) && .wi[1L] < length(args_raw)) {
  max(1L, as.integer(args_raw[.wi[1L] + 1L]))
} else {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || n < 1L) n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  max(1L, min(n %||% 4L, 8L))
}

# Everything that isn't a flag or the worker count is a family id
.flags <- c("--quick", "--parallel", "--workers",
            if (length(.wi)) args_raw[.wi[1L] + 1L] else character(0))
FAMILIES  <- args_raw[!args_raw %in% .flags]
if (length(FAMILIES) == 0L) FAMILIES <- c("OSSL_VisNIR", "OSSL_MIR")

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# ---- Configuration ----------------------------------------------------------

K_FOLDS   <- if (QUICK) 3L  else 10L
EPOCHS    <- if (QUICK) 20L else 80L
BATCH     <- 32L
PAT_ES    <- if (QUICK) 5L  else 10L
PAT_LR    <- if (QUICK) 3L  else 5L
LATENT    <- 16L
MIN_N     <- 50L
OUT_DIR   <- "models"
CACHE_DIR <- autoSpectra::ossl_cache_dir()

# In quick mode only train these 4 sentinel properties
QUICK_PROPS <- c("oc", "clay.tot", "ph.h2o", "n.tot")

if (QUICK)    message("*** QUICK MODE: K=", K_FOLDS, ", epochs=", EPOCHS,
                       ", properties=", paste(QUICK_PROPS, collapse = ", "))
if (PARALLEL) message("*** PARALLEL MODE: up to ", N_WORKERS, " concurrent workers")

# ---- Helper: LODO k-fold indices -------------------------------------------

#' @param dataset_col Character vector of dataset labels (one per sample)
#' @param k Number of folds
#' @param seed RNG seed
#' @return List of length k; each element is integer vector of test indices
stratified_kfold <- function(dataset_col, k = 10L, seed = 42L) {
  set.seed(seed)
  folds    <- vector("list", k)
  datasets <- unique(dataset_col)
  ds_order <- sample(datasets)
  ds_fold  <- setNames(((seq_along(ds_order) - 1L) %% k) + 1L, ds_order)
  for (fold in seq_len(k))
    folds[[fold]] <- which(ds_fold[dataset_col] == fold)

  empty <- which(vapply(folds, length, integer(1)) == 0L)
  if (length(empty) > 0L) {
    message("  Note: ", length(empty), " empty fold(s) — ",
            "using random assignment instead of dataset stratification.")
    n    <- length(dataset_col)
    idx  <- sample.int(n)
    cuts <- round(seq(0, n, length.out = k + 1L))
    for (fold in empty)
      folds[[fold]] <- idx[(cuts[fold] + 1L):cuts[fold + 1L]]
  }
  folds
}

# ---- Helper: single-property CV training ------------------------------------

train_property_cv <- function(X, y, dataset_col, family_id, prop,
                               k       = K_FOLDS,
                               epochs  = EPOCHS,
                               batch   = BATCH,
                               latent  = LATENT,
                               pat_es  = PAT_ES,
                               pat_lr  = PAT_LR,
                               out_dir = OUT_DIR) {
  n <- nrow(X)
  if (n < MIN_N) {
    message("    SKIP: need >= ", MIN_N, ", have ", n)
    return(NULL)
  }

  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  # K-fold LODO cross-validation
  folds      <- stratified_kfold(dataset_col, k = k)
  fold_rmse  <- numeric(k); fold_bias <- numeric(k)
  fold_r2    <- numeric(k); fold_rpiq <- numeric(k)
  fold_ccc   <- numeric(k)

  for (fold in seq_len(k)) {
    te_idx <- folds[[fold]]
    if (length(te_idx) == 0L) {
      fold_rmse[fold] <- NA; fold_bias[fold] <- NA; fold_r2[fold] <- NA
      fold_rpiq[fold] <- NA; fold_ccc[fold]  <- NA
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
    yhat_z <- .extract_prediction(mdl_fold, X_te)
    yhat   <- yhat_z * sd_y + mean_y
    mets   <- metrics_from_y(y_te, yhat)
    fold_rmse[fold] <- mets$RMSE; fold_bias[fold] <- mets$Bias
    fold_r2[fold]   <- mets$R2;   fold_rpiq[fold] <- mets$RPIQ
    fold_ccc[fold]  <- mets$CCC
    keras3::clear_session()
    rm(mdl_fold); gc(verbose = FALSE)
  }

  cv <- list(
    RMSE_mean = mean(fold_rmse, na.rm = TRUE),
    RMSE_sd   = stats::sd(fold_rmse, na.rm = TRUE),
    Bias_mean = mean(fold_bias, na.rm = TRUE),
    Bias_sd   = stats::sd(fold_bias, na.rm = TRUE),
    R2_mean   = mean(fold_r2,   na.rm = TRUE),
    R2_sd     = stats::sd(fold_r2,   na.rm = TRUE),
    RPIQ_mean = mean(fold_rpiq, na.rm = TRUE),
    RPIQ_sd   = stats::sd(fold_rpiq, na.rm = TRUE),
    CCC_mean  = mean(fold_ccc,  na.rm = TRUE),
    CCC_sd    = stats::sd(fold_ccc,  na.rm = TRUE),
    n_folds   = k,
    fold_RMSE = fold_rmse, fold_Bias = fold_bias,
    fold_R2   = fold_r2,   fold_RPIQ = fold_rpiq,
    fold_CCC  = fold_ccc
  )

  # --- Final model: train on ALL data, 80/10/10 split ----------------------
  sp      <- split_idx(n, seed = 42L, p_train = 0.8, p_cal = 0.1, p_test = 0.1)
  X_tr    <- X[sp$train, , drop = FALSE]; y_tr <- y_z[sp$train]
  X_ca    <- X[sp$calib, , drop = FALSE]; y_ca <- y[sp$calib]
  X_te    <- X[sp$test,  , drop = FALSE]; y_te <- y[sp$test]

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

  # Conformal calibration (absolute residuals on calib set)
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

  latent_stats <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = d_lat,
    thr95 = stats::qchisq(0.95, df = d_lat)
  )

  save_metrics_json(base,
                    metrics    = mets_te,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = feature_minmax(X),
                    latent     = latent_stats)

  # Save CV results for parallel aggregation
  saveRDS(list(cv = cv, test = mets_te, n = n),
          paste0(base, "_cv.rds"))

  list(cv = cv, test = mets_te, n = n)
}

# ---- Parallel dispatch (subprocess isolation for TF/Keras) ------------------
#
# Each property is trained in a separate R session via callr::r_bg() so that
# TensorFlow sessions do not interfere. Pre-processed X_proc and y values are
# serialised once by the master and read by each worker from a temp RDS file.

dispatch_parallel <- function(props, worker_data_rds, family_id,
                               max_workers, out_dir,
                               k, epochs, batch, pat_es, pat_lr, latent,
                               script_dir = getwd()) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    message("  'callr' not installed — falling back to sequential. ",
            "Install with: install.packages('callr')")
    return(NULL)
  }

  # The function each worker subprocess will run
  worker_fn <- function(prop, data_rds, script_dir,
                        k, epochs, batch, pat_es, pat_lr, latent,
                        out_dir, min_n) {
    setwd(script_dir)
    if (file.exists("DESCRIPTION")) {
      pkgload::load_all(".", quiet = TRUE)
    } else {
      library(autoSpectra)
    }
    suppressPackageStartupMessages(library(keras3))

    wdata      <- readRDS(data_rds)
    X_proc     <- wdata$X_proc
    ds_labels  <- wdata$dataset_labels
    fam_id     <- wdata$family_id
    y_df       <- wdata$y_df

    y    <- suppressWarnings(as.numeric(y_df[[prop]]))
    keep <- is.finite(y) & is.finite(rowSums(X_proc))
    n_ok <- sum(keep)

    message(sprintf("  [%-20s]  n=%-6d  training ...", prop, n_ok))

    if (n_ok < min_n) {
      message("    SKIP (n < ", min_n, ")")
      return(NULL)
    }

    train_property_cv(
      X           = X_proc[keep, , drop = FALSE],
      y           = y[keep],
      dataset_col = ds_labels[keep],
      family_id   = fam_id,
      prop        = prop,
      k           = k,      epochs  = epochs,
      batch       = batch,  latent  = latent,
      pat_es      = pat_es, pat_lr  = pat_lr,
      out_dir     = out_dir
    )
  }

  # --- Pool management: run max_workers at a time ---------------------------
  pool  <- list()   # name -> r_bg process
  done  <- list()   # name -> result
  queue <- as.list(props)

  message("  Dispatching ", length(props), " properties across ",
          max_workers, " parallel workers ...")

  while (length(queue) > 0L || length(pool) > 0L) {
    # Launch workers until pool is full
    while (length(queue) > 0L && length(pool) < max_workers) {
      prop  <- queue[[1L]]; queue <- queue[-1L]
      message("  -> launching worker: ", prop)
      pool[[prop]] <- callr::r_bg(
        func = worker_fn,
        args = list(
          prop       = prop,
          data_rds   = worker_data_rds,
          script_dir = script_dir,
          k          = k,      epochs  = epochs,
          batch      = batch,  pat_es  = pat_es,
          pat_lr     = pat_lr, latent  = latent,
          out_dir    = out_dir, min_n  = MIN_N
        ),
        package = FALSE
      )
    }

    # Harvest finished workers
    alive <- vapply(names(pool), function(nm) pool[[nm]]$is_alive(), logical(1))
    for (nm in names(pool)[!alive]) {
      result <- tryCatch(pool[[nm]]$get_result(), error = function(e) {
        message("  Worker ERROR [", nm, "]: ", conditionMessage(e))
        NULL
      })
      done[[nm]] <- result
      if (!is.null(result)) {
        cv <- result$cv
        message(sprintf(
          "  [%-20s] CV: R2=%.3f\u00b1%.3f  RMSE=%.3f\u00b1%.3f  Bias=%.3f\u00b1%.3f  CCC=%.3f\u00b1%.3f  RPIQ=%.2f\u00b1%.2f",
          nm, cv$R2_mean, cv$R2_sd, cv$RMSE_mean, cv$RMSE_sd,
          cv$Bias_mean, cv$Bias_sd, cv$CCC_mean, cv$CCC_sd,
          cv$RPIQ_mean, cv$RPIQ_sd))
      }
      pool[[nm]] <- NULL
    }
    pool <- pool[!vapply(names(pool), is.null, logical(1))]

    if (length(pool) > 0L) Sys.sleep(10)
  }

  done
}

# ---- Main training loop -----------------------------------------------------

for (fam_id in FAMILIES) {
  fam   <- get_family(fam_id)
  props <- if (QUICK) intersect(QUICK_PROPS, fam$properties) else fam$properties

  message("\n", strrep("=", 60))
  message("Training: ", fam$label)
  message(strrep("=", 60))
  message("Grid     : ", length(fam$wavegrid), " bands (",
          min(fam$wavegrid), " - ", max(fam$wavegrid), ")")
  message("Pipeline : ", paste(fam$preprocess, collapse = " -> "))
  message("CV folds : ", K_FOLDS, " (LODO — stratified by OSSL dataset)")
  message("Metrics  : RMSE | Bias | R\u00b2 | RPIQ | CCC")

  # ---- Load and prepare OSSL data ------------------------------------------
  message("\nLoading OSSL data ...")
  joined_df <- ossl_prepare(
    sensor_type         = fam$sensor_type,
    cache_dir           = CACHE_DIR,
    properties          = fam$properties,
    download_if_missing = TRUE
  )
  message("  Samples after join: ", nrow(joined_df))

  ds_col <- intersect(c("dataset.code_ascii_txt", "dataset_code", "location_id"),
                      names(joined_df))
  if (length(ds_col) == 0L) {
    message("  WARNING: no dataset column — using random CV folds")
    dataset_labels <- as.character(seq_len(nrow(joined_df)))
  } else {
    dataset_labels <- as.character(joined_df[[ds_col[1L]]])
    n_ds <- length(unique(dataset_labels))
    message("  Datasets found: ", n_ds,
            " (", paste(head(unique(dataset_labels), 5L), collapse = ", "),
            if (n_ds > 5L) ", ..." else "", ")")
  }

  # Resample to canonical grid
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)

  message("  Preprocessing ...")
  X_proc <- apply_pipeline(X_res, fam$preprocess)

  # ---- Parallel or sequential property training ----------------------------
  summary_list <- list()

  if (PARALLEL) {
    # Serialise preprocessed data once; workers read from this file
    tmp_rds <- tempfile("ossl_worker_", tmpdir = tempdir(), fileext = ".rds")
    saveRDS(list(
      X_proc         = X_proc,
      y_df           = joined_df[, c("Soil_ID", fam$properties), drop = FALSE],
      dataset_labels = dataset_labels,
      family_id      = fam_id
    ), tmp_rds)
    on.exit(unlink(tmp_rds), add = TRUE)

    # Filter to properties with enough samples
    props_to_run <- Filter(function(p) {
      y    <- suppressWarnings(as.numeric(joined_df[[p]]))
      keep <- is.finite(y) & is.finite(rowSums(X_proc))
      n_ok <- sum(keep)
      if (n_ok < MIN_N) {
        message(sprintf("  [%-20s]  SKIP  (n=%d)", p, n_ok))
        FALSE
      } else TRUE
    }, props)

    par_results <- dispatch_parallel(
      props            = props_to_run,
      worker_data_rds  = tmp_rds,
      family_id        = fam_id,
      max_workers      = N_WORKERS,
      out_dir          = OUT_DIR,
      k                = K_FOLDS,  epochs  = EPOCHS,
      batch            = BATCH,    pat_es  = PAT_ES,
      pat_lr           = PAT_LR,   latent  = LATENT,
      script_dir       = getwd()
    )

    if (is.null(par_results)) {
      PARALLEL <- FALSE   # callr unavailable — fall through to sequential below
    } else {
      for (prop in names(par_results)) {
        r <- par_results[[prop]]
        if (!is.null(r)) {
          cv <- r$cv
          summary_list[[prop]] <- list(
            property   = prop,
            label      = property_label(prop),
            n          = r$n,
            R2_mean    = cv$R2_mean,   R2_sd    = cv$R2_sd,
            RMSE_mean  = cv$RMSE_mean, RMSE_sd  = cv$RMSE_sd,
            Bias_mean  = cv$Bias_mean, Bias_sd  = cv$Bias_sd,
            RPIQ_mean  = cv$RPIQ_mean, RPIQ_sd  = cv$RPIQ_sd,
            CCC_mean   = cv$CCC_mean,  CCC_sd   = cv$CCC_sd
          )
        }
      }
    }
  }

  # Sequential fallback (or when --parallel not requested)
  if (!PARALLEL) {
    for (prop in props) {
      y    <- suppressWarnings(as.numeric(joined_df[[prop]]))
      keep <- is.finite(y) & is.finite(rowSums(X_proc))
      n_ok <- sum(keep)

      if (n_ok < MIN_N) {
        message(sprintf("  [%-20s]  SKIP  (n=%d)", prop, n_ok))
        next
      }

      message(sprintf("  [%-20s]  n=%-6d  training ...", prop, n_ok))

      result <- tryCatch(
        train_property_cv(
          X           = X_proc[keep, , drop = FALSE],
          y           = y[keep],
          dataset_col = dataset_labels[keep],
          family_id   = fam_id,
          prop        = prop
        ),
        error = function(e) {
          message("    ERROR: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(result)) {
        cv <- result$cv
        message(sprintf(
          "    CV: R2=%.3f\u00b1%.3f  RMSE=%.3f\u00b1%.3f  Bias=%.3f\u00b1%.3f  CCC=%.3f\u00b1%.3f  RPIQ=%.2f\u00b1%.2f",
          cv$R2_mean, cv$R2_sd, cv$RMSE_mean, cv$RMSE_sd,
          cv$Bias_mean, cv$Bias_sd, cv$CCC_mean, cv$CCC_sd,
          cv$RPIQ_mean, cv$RPIQ_sd))
        summary_list[[prop]] <- list(
          property   = prop,
          label      = property_label(prop),
          n          = result$n,
          R2_mean    = cv$R2_mean,   R2_sd    = cv$R2_sd,
          RMSE_mean  = cv$RMSE_mean, RMSE_sd  = cv$RMSE_sd,
          Bias_mean  = cv$Bias_mean, Bias_sd  = cv$Bias_sd,
          RPIQ_mean  = cv$RPIQ_mean, RPIQ_sd  = cv$RPIQ_sd,
          CCC_mean   = cv$CCC_mean,  CCC_sd   = cv$CCC_sd
        )
      }

      keras3::clear_session()
      gc(verbose = FALSE)
    }
  }

  # ---- Save family-level metrics summary -----------------------------------
  summary_path <- file.path(OUT_DIR, fam_id, "metrics_summary.json")
  dir_create(dirname(summary_path))
  jsonlite::write_json(list(
    family       = fam_id,
    ossl_version = fam$ossl_version,
    k_folds      = K_FOLDS,
    strategy     = "LODO — stratified by OSSL contributing dataset",
    metrics      = c("RMSE", "Bias", "R2", "RPIQ", "CCC"),
    references   = list(
      RPIQ = "Bellon-Maurel et al. (2010)",
      CCC  = "Lin (1989)",
      LODO = "Safanelli et al. (2023)"
    ),
    properties   = summary_list
  ), summary_path, auto_unbox = TRUE, pretty = TRUE)

  message("\nSummary saved: ", summary_path)
  message("Done: ", fam_id, " (", length(summary_list), " properties trained)\n")
}

message("\n", strrep("=", 60))
message("All families complete. Run run_autoSpectra() to launch the app.")
message(strrep("=", 60))
