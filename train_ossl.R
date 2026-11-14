# train_ossl.R - Training orchestrator for autoSpectra (Safanelli 2025 faithful)
#
# Trains four OSSL model families:
#   OSSL_VisNIR        - masked soilVAE on 350-2500 nm
#   OSSL_MIR           - masked soilVAE on 600-4000 cm-1
#   OSSL_VisNIR_Cubist - Cubist + PCA(120) + SNV (Safanelli replication)
#   OSSL_MIR_Cubist    - Cubist + PCA(120) + SNV (Safanelli replication)
#   OSSL_VisNIR_PCAVAE - SNV + PCA(120) + soilVAE ensemble (only the model differs from Cubist)
#   OSSL_MIR_PCAVAE    - SNV + PCA(120) + soilVAE ensemble
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
#   models/<family_id>/cv10_summary.json  (random 10-fold CV, --cv10)
#   models/<family_id>/lodo_summary.json  (leave-one-dataset-out, --cv10 --lodo)
#
# Usage:
#   Rscript train_ossl.R                            # all 4 families, full
#   Rscript train_ossl.R OSSL_MIR                   # one family only
#   Rscript train_ossl.R --only-vae                 # both VAE families
#   Rscript train_ossl.R --only-cubist              # both Cubist families
#   Rscript train_ossl.R --only-pcavae              # both PCA-VAE families
#   Rscript train_ossl.R --cv10                     # 10-fold CV refit mode
#   Rscript train_ossl.R --quick                    # 3-fold, 20 epochs, 4 props
#   Rscript train_ossl.R --props oc,clay.tot,ph.h2o # selected properties
#   Rscript train_ossl.R OSSL_NeoSpectra_Cubist     # OSSL NeoSpectra library (no African soils)
#   Rscript train_ossl.R --subset kssl --out-dir models_kssl  # KSSL rows only (Safanelli's "kssl" models)
#   Rscript train_ossl.R --cv10 --pca-fit library --pca-scale  # PCA on every spectrum
#       of the library except the test fold, bands scaled (OSSL-style PCA, no leak)

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
FORCE       <- has_flag("--force")   # retrain properties that already finished
LODO        <- has_flag("--lodo")
ONLY_VAE    <- has_flag("--only-vae")
ONLY_CUBIST <- has_flag("--only-cubist")
ONLY_PCAVAE <- has_flag("--only-pcavae")

SUBSET    <- arg_val("--subset", "ossl")      # "ossl" (all rows) or "kssl"
if (!SUBSET %in% c("ossl", "kssl")) stop("--subset must be ossl or kssl")
# Rows the PCA of the Cubist / PCA-VAE families is fitted on:
# "train" = the training rows that have the property (default);
# "library" = every spectrum of the sensor except the held-out samples (the
# test fold in --cv10, the 20% test split otherwise).
PCA_FIT   <- arg_val("--pca-fit", "train")
if (!PCA_FIT %in% c("train", "library")) stop("--pca-fit must be train or library")
PCA_SCALE <- has_flag("--pca-scale")
PROPS_ARG <- arg_val("--props", NULL)
PROPS     <- if (!is.null(PROPS_ARG)) trimws(strsplit(PROPS_ARG, ",")[[1]]) else NULL

# --workers for parallelism (currently unused: keras + cubist are sequential
# per-property to avoid GPU contention; retained for future parallel support)
N_WORKERS <- suppressWarnings(as.integer(arg_val("--workers", "1")))
if (is.na(N_WORKERS)) N_WORKERS <- 1L

# Resolve family list from positional args (everything not starting with --
# and not the value of a known flag)
flag_values <- c(arg_val("--props"), arg_val("--workers"), arg_val("--subset"),
                 arg_val("--out-dir"), arg_val("--pca-fit"))
ALL_FAMILIES <- c("OSSL_VisNIR", "OSSL_MIR",
                  "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist",
                  "OSSL_VisNIR_PCAVAE", "OSSL_MIR_PCAVAE")
POS <- args_raw[!grepl("^--", args_raw) & !(args_raw %in% flag_values)]
FAMILIES <- if (length(POS) > 0L) POS else ALL_FAMILIES
if (ONLY_VAE)    FAMILIES <- intersect(FAMILIES, c("OSSL_VisNIR", "OSSL_MIR"))
if (ONLY_CUBIST) FAMILIES <- intersect(FAMILIES,
                                       c("OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist"))
if (ONLY_PCAVAE) FAMILIES <- intersect(FAMILIES,
                                       c("OSSL_VisNIR_PCAVAE", "OSSL_MIR_PCAVAE"))

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
OUT_DIR   <- arg_val("--out-dir", "models")
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
message("  PCA fit       : ", PCA_FIT, if (PCA_SCALE) ", bands scaled" else "")

# ---- Per-family training -------------------------------------------------

.safe_nm <- function(family_id, mode = c("single", "cv10")) {
  mode <- match.arg(mode)
  dir_create(file.path(OUT_DIR, family_id))
  # LODO and random 10-fold CV answer different questions (unseen datasets vs
  # Safanelli's S2 protocol), so they must not overwrite each other.
  file.path(OUT_DIR, family_id,
            if (mode == "single") "metrics_summary.json"
            else if (LODO) "lodo_summary.json" else "cv10_summary.json")
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
  jsonlite::write_json(summary, out_path, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA)
  message("  Summary saved: ", out_path)
}

# Keras/TF accumulates graph and thread-pool state across successive models in
# one process. On macOS that surfaced as a hard C++ abort when building the
# SECOND model of a run, after the first had trained and saved cleanly:
#
#   libc++abi: terminating due to uncaught exception of type
#   std::__1::system_error: condition_variable wait failed: Invalid argument
#
# clear_session() drops the backend state between properties; gc() then releases
# the R-side references to the freed Python objects.
.reset_keras_session <- function(algo) {
  if (identical(algo, "cubist_pca")) return(invisible(NULL))
  if (!requireNamespace("keras3", quietly = TRUE)) return(invisible(NULL))
  try(keras3::clear_session(), silent = TRUE)
  invisible(gc(verbose = FALSE))
}

# A property counts as finished only when BOTH its model and its metrics are on
# disk; a model without metrics means the run died mid-save and must be redone.
.property_done <- function(family_id, prop, algo) {
  base <- file.path(OUT_DIR, family_id, "models", prop)
  model_ok <- if (identical(algo, "cubist_pca"))
    file.exists(paste0(base, "_cubist.rds"))
  else if (identical(algo, "pca_vae"))
    file.exists(paste0(base, "_pcavae.rds"))
  else
    file.exists(paste0(base, ".keras")) || file.exists(paste0(base, ".h5"))
  model_ok && file.exists(paste0(base, "_metrics.json"))
}

# Training inputs come from the package (.family_training_inputs in
# R/training_inputs.R), which uses the same preprocessing call as
# predict_soil(). Targets go through apply_target_transform(), which also drops
# physically impossible lab values.

# Library for the PCA of the Cubist / PCA-VAE families (--pca-fit library):
# every spectrum of the sensor, labelled or not. NULLs otherwise.
.pca_library <- function(fam, pp) {
  if (PCA_FIT != "library" || !fam$algorithm %in% c("cubist_pca", "pca_vae"))
    return(list(X = NULL, groups = NULL))
  message("  PCA library: ", sum(pp$has_spectra), " spectra",
          if (PCA_SCALE) ", bands scaled" else "")
  list(X = pp$X[pp$has_spectra, , drop = FALSE], groups = pp$groups[pp$has_spectra])
}

.train_family_single <- function(family_id, props, joined_df) {
  fam  <- get_family(family_id)
  algo <- fam$algorithm %||% "vae"

  pp   <- .family_training_inputs(fam, joined_df)
  lib  <- .pca_library(fam, pp)

  results <- list()
  for (prop in props) {
    # Resume: a finished property leaves both its model and its metrics behind.
    # Skipping it lets a crashed run be restarted without redoing hours of work
    # (see the TF condition_variable crash noted in .reset_keras_session()).
    if (!FORCE && .property_done(family_id, prop, algo)) {
      message("  Done already [", prop, "] - skipping (use --force to retrain)")
      results[[prop]] <- TRUE
      next
    }

    y    <- apply_target_transform(joined_df[[prop]], prop)
    keep <- is.finite(y) & pp$has_spectra
    if (sum(keep) < MIN_N) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      next
    }
    message("  Training [", prop, "] (n=", sum(keep), ", sample groups=",
            length(unique(pp$groups[keep])), ") ...")
    res <- try({
      if (identical(algo, "vae_masked")) {
        train_maskedVAE(pp$X_raw[keep, , drop = FALSE], y[keep], family_id, prop,
                        mask        = pp$mask[keep, , drop = FALSE],
                        groups      = pp$groups[keep],
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
        train_cubist_pca(pp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                         out_dir = OUT_DIR, min_n = MIN_N,
                         groups = pp$groups[keep],
                         X_lib = lib$X, lib_groups = lib$groups, pca_scale = PCA_SCALE)
      } else if (identical(algo, "pca_vae")) {
        train_pca_vae(pp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                      out_dir = OUT_DIR, min_n = MIN_N,
                      groups = pp$groups[keep],
                      X_lib = lib$X, lib_groups = lib$groups, pca_scale = PCA_SCALE,
                      epochs = if (QUICK) 20L else 200L)
      } else {
        train_soilVAE(pp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                      out_dir = OUT_DIR, groups = pp$groups[keep],
                      epochs = EPOCHS, batch_size = BATCH,
                      patience_es = PAT_ES, patience_rl = PAT_LR,
                      latent_dim = LATENT, kl_beta = KL_BETA,
                      kl_warmup_epochs = KL_WARMUP, min_n = MIN_N)
      }
    }, silent = FALSE)
    results[[prop]] <- res
    .reset_keras_session(algo)
  }
  results
}

.train_family_cv10 <- function(family_id, props, joined_df) {
  fam  <- get_family(family_id)
  pp   <- .family_training_inputs(fam, joined_df)
  masked <- identical(fam$algorithm, "vae_masked")
  X    <- if (masked) pp$X_raw else pp$X
  lib   <- .pca_library(fam, pp)
  X_lib <- lib$X; g_lib <- lib$groups

  cv_summaries <- list()
  # One file per property, so a run of several days can be resumed after a
  # crash; the family summary is assembled from these files.
  prop_dir <- file.path(OUT_DIR, family_id, if (LODO) "lodo" else "cv10")
  dir_create(prop_dir)
  for (prop in props) {
    pf <- file.path(prop_dir, paste0(prop, ".json"))
    if (!FORCE && file.exists(pf)) {
      message("  Done already [", prop, "] - ", pf)
      cv_summaries[[prop]] <- jsonlite::read_json(pf, simplifyVector = FALSE)
      next
    }
    y    <- apply_target_transform(joined_df[[prop]], prop)
    keep <- is.finite(y) & pp$has_spectra
    if (sum(keep) < MIN_N) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      next
    }
    y_k  <- y[keep]; X_k <- X[keep, , drop = FALSE]
    g_k  <- pp$groups[keep]

    folds <- if (LODO &&
                 "dataset.code_ascii_txt" %in% names(joined_df)) {
      ds <- joined_df[["dataset.code_ascii_txt"]][keep]
      stratified_kfold(ds, k = K_FOLDS, seed = 42L)
    } else {
      random_kfold(length(y_k), k = K_FOLDS, seed = 42L, groups = g_k)
    }

    message("  10-fold CV [", prop, "] (n=", length(y_k), ") ...")
    cv_res <- try(
      cv_train_property(X_k, y_k, family_id, prop, folds, fam,
                        mask   = if (masked) pp$mask[keep, , drop = FALSE],
                        groups = g_k,
                        X_lib = X_lib, lib_groups = g_lib, pca_scale = PCA_SCALE,
                        epochs = EPOCHS, batch_size = BATCH,
                        patience_es = PAT_ES, patience_rl = PAT_LR,
                        latent_dim = LATENT, kl_beta = KL_BETA,
                        kl_warmup_epochs = KL_WARMUP),
      silent = FALSE
    )
    cv_summaries[[prop]] <- if (inherits(cv_res, "try-error")) NULL else cv_res
    if (!is.null(cv_summaries[[prop]]))
      jsonlite::write_json(cv_summaries[[prop]], pf, auto_unbox = TRUE, pretty = TRUE,
                           digits = NA)
    .reset_keras_session(fam$algorithm %||% "vae")
    ms <- cv_summaries[[prop]]$mean_sd
    if (!is.null(ms)) {
      message(sprintf(
        "    RMSE=%.3f+-%.3f  R2=%.3f+-%.3f  RPIQ=%.2f+-%.2f  CCC=%.3f+-%.3f",
        ms$RMSE_mean, ms$RMSE_sd, ms$R2_mean, ms$R2_sd,
        ms$RPIQ_mean, ms$RPIQ_sd, ms$CCC_mean, ms$CCC_sd))
      po <- cv_summaries[[prop]]$pooled
      if (!is.null(po))
        message(sprintf("    pooled: R2=%.3f  cor2=%.3f  RMSE=%.3f  (n=%d)",
                        po$R2, po$cor2, po$RMSE, po$n))
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
    if (SUBSET == "kssl") {
      if (sensor_type == "neospectra") stop("--subset kssl applies to visnir and mir only")
      df <- df[df[["dataset.code_ascii_txt"]] %in% "KSSL.SSL", , drop = FALSE]
      message("  KSSL subset: ", nrow(df), " rows")
    }
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

  props <- intersect(target_props, fam$properties)

  # Skip a fully-trained family BEFORE touching the data. Loading and
  # preprocessing a sensor costs minutes (132k x 1076 for VisNIR, x1701 for
  # MIR), and the supervisor relaunches this script after every TF abort, so
  # deciding first and loading second is what makes a restart cheap.
  if (!CV10 && !FORCE && length(props) > 0L &&
      all(vapply(props, function(p)
            .property_done(fid, p, fam$algorithm %||% "vae"), logical(1)))) {
    message("  All ", length(props),
            " requested properties already trained - skipping family")
    .summarise_family(fid, setNames(as.list(rep(TRUE, length(props))), props))
    next
  }

  joined_df <- .load_joined(fam$sensor_type)

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
