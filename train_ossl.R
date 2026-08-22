# train_ossl.R — autoSpectra OSSL model training script
#
# Trains one soilVAE model per soil property × spectral domain using the full
# OSSL v1.2 corpus. Two sensor-agnostic families:
#
#   OSSL_VisNIR — Masked VAE: encoder sees [SNV-spectrum | binary-mask].
#                 Trained with 3 spectral coverage profiles (full, NIR, SW-NIR),
#                 so a single model handles any VisNIR sensor.
#
#   OSSL_MIR    — Standard VAE: all MIR instruments share the 600-4000 cm⁻¹ grid.
#
# Validation: Grouped repeated stratified K-fold by Soil_ID.
#   Replicate spectra stay together; target quantiles are balanced across folds.
#
# Metrics: RMSE | Bias | R² | RPIQ (Bellon-Maurel 2010) | CCC (Lin 1989)
#
# Outputs (per family):
#   models/<family_id>/models/<prop>.keras        — soilVAE weights
#   models/<family_id>/models/<prop>_scaler.rds   — z-score parameters
#   models/<family_id>/models/<prop>_cv.rds       — CV fold results
#   models/<family_id>/metrics_summary.json       — all properties, mean±SD
#
# Usage:
#   Rscript train_ossl.R                          # both families, all properties
#   Rscript train_ossl.R OSSL_VisNIR              # VisNIR only
#   Rscript train_ossl.R OSSL_MIR                 # MIR only
#   Rscript train_ossl.R --props oc,clay.tot      # specific properties
#   Rscript train_ossl.R --quick                  # 3-fold, 20 epochs, oc only
#   Rscript train_ossl.R --sequential             # disable parallel
#   Rscript train_ossl.R --workers 4              # max parallel workers
#   Rscript train_ossl.R --seed 42                # reproducible results
#   Rscript train_ossl.R --epochs 100             # override epoch cap
#   Rscript train_ossl.R --sample-size 2000       # limit training set size

# ---- Bootstrap ---------------------------------------------------------------

load_training_context <- function(script_dir = getwd()) {
  python_has_keras <- function(python) {
    if (is.null(python) || !nzchar(python) || !file.exists(python))
      return(FALSE)
    code <- paste(
      "import importlib.util, sys",
      "mods = ('keras','tensorflow','numpy')",
      "missing=[m for m in mods if importlib.util.find_spec(m) is None]",
      "sys.exit(1 if missing else 0)",
      sep = "\n"
    )
    out <- tryCatch(
      suppressWarnings(system2(python, c("-c", shQuote(code)),
                               stdout = TRUE, stderr = TRUE)),
      error = function(e) structure(character(0), status = 1L)
    )
    status <- attr(out, "status")
    if (is.null(status) || length(status) == 0L) status <- 0L
    identical(status, 0L)
  }

  bootstrap_python <- function() {
    env_py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
    if (python_has_keras(env_py)) return(invisible(env_py))

    home   <- Sys.getenv("USERPROFILE",
                         unset = Sys.getenv("HOME", unset = ""))
    venv_roots <- unique(c(
      file.path(home, "OneDrive - University of Florida",
                "Documentos", ".virtualenvs"),
      file.path(home, "OneDrive - University of Florida",
                "Documents", ".virtualenvs"),
      file.path(home, ".virtualenvs")
    ))
    candidates <- character(0)
    for (root in venv_roots) {
      for (venv in c("r-tf311", "r-tensorflow", "r-tf", "r-reticulate")) {
        candidates <- c(candidates,
                        file.path(root, venv, "Scripts", "python.exe"),
                        file.path(root, venv, "bin", "python"))
      }
    }
    py_bin <- Sys.which("python3")
    if (!nzchar(py_bin)) py_bin <- Sys.which("python")
    if (nzchar(py_bin)) candidates <- c(candidates, py_bin)
    candidates <- unique(candidates[nzchar(candidates)])

    for (py in candidates) {
      if (!python_has_keras(py)) next
      Sys.setenv(RETICULATE_PYTHON = py)
      return(invisible(py))
    }
    invisible(NULL)
  }

  bootstrap_python()

  desc_path <- file.path(script_dir, "DESCRIPTION")
  if (file.exists(desc_path)) {
    r_files <- sort(list.files(
      file.path(script_dir, "R"),
      pattern = "\\.[Rr]$",
      full.names = TRUE
    ))
    for (path in r_files) sys.source(path, envir = globalenv())
  } else {
    suppressPackageStartupMessages(library(autoSpectra))
  }
  invisible(TRUE)
}

message("train_ossl.R: loading model code ...")
load_training_context(getwd())
message("train_ossl.R: model code loaded")

# ---- Argument parsing --------------------------------------------------------

if (!exists("args_raw", inherits = FALSE))
  args_raw <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

QUICK      <- "--quick"       %in% args_raw
PARALLEL   <- !"--sequential" %in% args_raw

# --workers N
.wi <- which(args_raw == "--workers")
N_WORKERS <- if (length(.wi) && .wi[1L] < length(args_raw)) {
  max(1L, as.integer(args_raw[.wi[1L] + 1L]))
} else {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || n < 1L)
    n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  max(1L, (n %||% 4L) - 2L)
}

# --props oc,clay.tot
.pi <- which(args_raw == "--props")
PROPS_FILTER <- if (length(.pi) && .pi[1L] < length(args_raw)) {
  strsplit(args_raw[.pi[1L] + 1L], ",", fixed = TRUE)[[1L]]
} else NULL

# --sample-size N
.ni <- which(args_raw == "--sample-size")
SAMPLE_N <- if (length(.ni) && .ni[1L] < length(args_raw)) {
  as.integer(args_raw[.ni[1L] + 1L])
} else NULL

# --seed N
.si <- which(args_raw == "--seed")
SEED <- if (length(.si) && .si[1L] < length(args_raw)) {
  as.integer(args_raw[.si[1L] + 1L])
} else NULL

# --epochs N
.ei <- which(args_raw == "--epochs")
EPOCHS_OVERRIDE <- if (length(.ei) && .ei[1L] < length(args_raw)) {
  as.integer(args_raw[.ei[1L] + 1L])
} else NULL

# Positional arguments = family IDs
.known_flags <- c(
  "--quick", "--sequential",
  "--workers",   if (length(.wi)) args_raw[.wi[1L] + 1L] else character(0),
  "--props",     if (length(.pi)) args_raw[.pi[1L] + 1L] else character(0),
  "--sample-size", if (length(.ni)) args_raw[.ni[1L] + 1L] else character(0),
  "--seed",      if (length(.si)) args_raw[.si[1L] + 1L] else character(0),
  "--epochs",    if (length(.ei)) args_raw[.ei[1L] + 1L] else character(0)
)
FAMILIES <- args_raw[!args_raw %in% .known_flags]
if (length(FAMILIES) == 0L) FAMILIES <- c("OSSL_VisNIR", "OSSL_MIR")

# ---- Training configuration --------------------------------------------------

K_FOLDS    <- if (QUICK) 3L  else 10L
EPOCHS     <- EPOCHS_OVERRIDE %||% if (QUICK) 20L else 80L
BATCH      <- 32L
PAT_ES     <- if (QUICK) 5L  else 10L
PAT_LR     <- if (QUICK) 3L  else 5L
LATENT     <- 16L
KL_BETA    <- 1e-3            # meaningful VAE regularisation (4e-5 is too small)
KL_WARMUP  <- if (QUICK) 10L else 30L
CV_REPEATS <- if (QUICK) 3L  else 1L
MIN_N      <- 50L
OUT_DIR    <- "models"
CACHE_DIR  <- ossl_cache_dir()
RUN_SEED   <- SEED %||% if (QUICK) 42L else NULL
QUICK_PROPS <- "oc"

message("  Families : ", paste(FAMILIES, collapse = ", "))
message("  K-folds  : ", K_FOLDS,
        "  epochs=", EPOCHS,
        "  latent=", LATENT,
        "  kl_beta=", KL_BETA)
if (QUICK)    message("*** QUICK MODE — sentinel: ", QUICK_PROPS)
if (!PARALLEL) message("*** SEQUENTIAL MODE")
if (!is.null(PROPS_FILTER))
  message("*** PROPS FILTER: ", paste(PROPS_FILTER, collapse = ", "))
if (!is.null(SAMPLE_N))
  message("*** SAMPLE SIZE: ", SAMPLE_N, " unique samples")
if (!is.null(RUN_SEED))
  message("*** SEED: ", RUN_SEED)
if (PARALLEL)
  message("*** PARALLEL: up to ", N_WORKERS, " workers")

# ---- Helper: format CV summary line -----------------------------------------

format_cv_summary <- function(cv, prefix = "CV") {
  sprintf(
    "%s: R2=%.3f\u00b1%.3f  RMSE=%.3f\u00b1%.3f  Bias=%.3f\u00b1%.3f  CCC=%.3f\u00b1%.3f  RPIQ=%.2f\u00b1%.2f",
    prefix,
    cv$R2_mean,   cv$R2_sd,
    cv$RMSE_mean, cv$RMSE_sd,
    cv$Bias_mean, cv$Bias_sd,
    cv$CCC_mean,  cv$CCC_sd,
    cv$RPIQ_mean, cv$RPIQ_sd
  )
}

format_cv_fold_lines <- function(cv, indent = "      ") {
  fs <- cv$fold_summary
  if (is.null(fs) || nrow(fs) < 1L) return(character(0))
  vapply(seq_len(nrow(fs)), function(i) {
    sprintf(
      "%sfold %d [%s, n=%d]: R2=%.3f  RMSE=%.3f  Bias=%.3f  CCC=%.3f  RPIQ=%.2f",
      indent, fs$fold[i], fs$held_out_dataset[i], fs$n_test[i],
      fs$R2[i], fs$RMSE[i], fs$Bias[i], fs$CCC[i], fs$RPIQ[i]
    )
  }, character(1))
}

# ---- Helper: parallel dispatch ----------------------------------------------

dispatch_parallel <- function(props, worker_data_rds, family_id,
                               max_workers, out_dir,
                               k, epochs, batch, pat_es, pat_lr,
                               latent, kl_beta, kl_warmup, cv_repeats,
                               seed, script_dir = getwd()) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    message("  'callr' not installed — using sequential mode. ",
            "Install with: install.packages('callr')")
    return(NULL)
  }

  worker_fn <- function(prop, data_rds, script_dir,
                        k, epochs, batch, pat_es, pat_lr,
                        latent, kl_beta, kl_warmup, cv_repeats,
                        seed, out_dir, min_n) {
    setwd(script_dir)

    # Minimal Python bootstrap
    python_has_keras <- function(python) {
      if (is.null(python) || !nzchar(python) || !file.exists(python))
        return(FALSE)
      code <- paste(
        "import importlib.util,sys",
        "mods=('keras','tensorflow','numpy')",
        "sys.exit(1 if [m for m in mods if importlib.util.find_spec(m) is None] else 0)",
        sep = "\n"
      )
      out <- tryCatch(
        suppressWarnings(system2(python, c("-c", shQuote(code)),
                                 stdout = TRUE, stderr = TRUE)),
        error = function(e) structure(character(0), status = 1L)
      )
      status <- attr(out, "status")
      if (is.null(status) || length(status) == 0L) status <- 0L
      identical(status, 0L)
    }
    env_py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
    if (!python_has_keras(env_py)) {
      home <- Sys.getenv("USERPROFILE",
                         unset = Sys.getenv("HOME", unset = ""))
      for (root in file.path(home, c(
        "OneDrive - University of Florida/Documentos/.virtualenvs",
        "OneDrive - University of Florida/Documents/.virtualenvs",
        ".virtualenvs"
      ))) {
        for (venv in c("r-tf311", "r-tensorflow", "r-tf", "r-reticulate")) {
          for (suffix in c("Scripts/python.exe", "bin/python")) {
            py <- file.path(root, venv, suffix)
            if (python_has_keras(py)) {
              Sys.setenv(RETICULATE_PYTHON = py)
              break
            }
          }
        }
      }
    }

    desc_path <- file.path(script_dir, "DESCRIPTION")
    if (file.exists(desc_path)) {
      r_files <- sort(list.files(file.path(script_dir, "R"),
                                 pattern = "\\.[Rr]$", full.names = TRUE))
      for (path in r_files) sys.source(path, envir = globalenv())
    } else {
      suppressPackageStartupMessages(library(autoSpectra))
    }

    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

    wdata         <- readRDS(data_rds)
    X_proc        <- wdata$X_proc
    ds_labels     <- wdata$dataset_labels
    sample_ids    <- wdata$sample_ids
    fam           <- wdata$family_config
    y_df          <- wdata$y_df
    mask_aug      <- wdata$mask
    n_sensors     <- wdata$n_sensors %||% 1L
    cv_eval_mask  <- wdata$cv_eval_mask

    y <- suppressWarnings(as.numeric(y_df[[prop]]))
    if (n_sensors > 1L) y <- rep(y, times = n_sensors)
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
      sample_id   = sample_ids[keep],
      cv_eval_mask = if (!is.null(cv_eval_mask)) cv_eval_mask[keep] else NULL,
      family_id   = fam$id,
      prop        = prop,
      k           = k,        epochs    = epochs,
      batch       = batch,    latent    = latent,
      pat_es      = pat_es,   pat_lr    = pat_lr,
      kl_beta     = kl_beta,  kl_warmup = kl_warmup,
      cv_repeats  = cv_repeats,
      seed        = seed,
      out_dir     = out_dir,
      min_n       = min_n,
      mask        = if (!is.null(mask_aug)) mask_aug[keep, , drop = FALSE] else NULL
    )
  }

  pool  <- list()
  done  <- list()
  queue <- as.list(props)

  message("  Dispatching ", length(props), " properties across ",
          max_workers, " parallel workers ...")

  while (length(queue) > 0L || length(pool) > 0L) {
    while (length(queue) > 0L && length(pool) < max_workers) {
      prop  <- queue[[1L]]; queue <- queue[-1L]
      message("  -> launching worker: ", prop)
      pool[[prop]] <- callr::r_bg(
        func = worker_fn,
        args = list(
          prop       = prop,
          data_rds   = worker_data_rds,
          script_dir = script_dir,
          k          = k,        epochs    = epochs,
          batch      = batch,    pat_es    = pat_es,
          pat_lr     = pat_lr,   latent    = latent,
          kl_beta    = kl_beta,  kl_warmup = kl_warmup,
          cv_repeats = cv_repeats, seed = seed,
          out_dir    = out_dir,  min_n = MIN_N
        ),
        package = FALSE
      )
    }

    alive <- vapply(names(pool), function(nm) pool[[nm]]$is_alive(), logical(1))
    for (nm in names(pool)[!alive]) {
      result <- tryCatch(
        pool[[nm]]$get_result(),
        error = function(e) {
          message("  Worker ERROR [", nm, "]: ", conditionMessage(e))
          NULL
        }
      )
      done[[nm]] <- result
      if (!is.null(result)) {
        cv <- result$cv
        message(sprintf("  [%-20s] %s", nm, format_cv_summary(cv)))
        for (line in format_cv_fold_lines(cv)) message(line)
      }
      pool[[nm]] <- NULL
    }
    pool <- pool[!vapply(names(pool), is.null, logical(1))]
    if (length(pool) > 0L) Sys.sleep(10)
  }

  done
}

# ---- Main training loop ------------------------------------------------------

for (fam_id in FAMILIES) {
  fam <- get_family(fam_id)

  props <- if (!is.null(PROPS_FILTER)) {
    intersect(PROPS_FILTER, fam$properties)
  } else if (QUICK) {
    intersect(QUICK_PROPS, fam$properties)
  } else {
    fam$properties
  }

  message("\n", strrep("=", 60))
  message("Training: ", fam$label)
  pipeline_label <- paste(fam$preprocess, collapse = " -> ")
  if (isTRUE(fam$masked))
    pipeline_label <- paste(pipeline_label, "SNV(masked)", sep = " -> ")
  message("Grid    : ", length(fam$wavegrid), " bands (",
          min(fam$wavegrid), " - ", max(fam$wavegrid), ")")
  message("Pipeline: ", pipeline_label)
  message("Metrics : RMSE | Bias | R\u00b2 | RPIQ | CCC")
  message(strrep("=", 60))

  # ---- Load OSSL data -------------------------------------------------------
  message("\nLoading OSSL data ...")
  joined_df <- ossl_prepare(
    sensor_type         = fam$sensor_type,
    cache_dir           = CACHE_DIR,
    properties          = fam$properties,
    download_if_missing = TRUE
  )
  message("  Samples after join: ", nrow(joined_df))

  ds_col <- intersect(c("dataset.code_ascii_txt", "dataset_code"), names(joined_df))
  dataset_labels <- if (length(ds_col) > 0L) {
    n_ds <- length(unique(joined_df[[ds_col[1L]]]))
    message("  Datasets: ", n_ds)
    as.character(joined_df[[ds_col[1L]]])
  } else {
    message("  WARNING: no dataset column — using single library label")
    rep("OSSL", nrow(joined_df))
  }

  sample_ids <- if ("Soil_ID" %in% names(joined_df)) {
    as.character(joined_df[["Soil_ID"]])
  } else {
    message("  WARNING: no Soil_ID column — using row ids as CV groups")
    as.character(seq_len(nrow(joined_df)))
  }

  # Optional sample size limit (replicate-aware)
  if (!is.null(SAMPLE_N)) {
    subset_info <- select_ossl_training_subset(
      joined_df,
      sample_n      = SAMPLE_N,
      sample_id_col = "Soil_ID",
      dataset_col   = if (length(ds_col)) ds_col[1L] else NULL,
      seed          = RUN_SEED %||% 42L
    )
    joined_df      <- subset_info$data
    dataset_labels <- dataset_labels[match(
      as.character(joined_df[["Soil_ID"]]),
      as.character(joined_df[["Soil_ID"]])
    )]
    dataset_labels <- as.character(joined_df[[ds_col[1L]]])
    sample_ids     <- as.character(joined_df[["Soil_ID"]])
    message("  Sample subset: ", subset_info$n_selected_rows, " rows")
  }

  # Limit k-folds to number of groups
  n_groups <- length(unique(sample_ids))
  ACTUAL_K <- max(1L, min(K_FOLDS, n_groups))
  if (ACTUAL_K != K_FOLDS)
    message("  Fold override: ", ACTUAL_K,
            " folds (limited by unique Soil_ID groups)")
  message("  CV folds: ", ACTUAL_K,
          " (grouped stratified K-fold by Soil_ID)  repeats=", CV_REPEATS)

  # ---- Preprocess ----------------------------------------------------------
  message("  Preprocessing ...")
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)
  X_proc  <- apply_pipeline(X_res, fam$preprocess)

  # ---- Spectral dropout augmentation for VisNIR masked family ---------------
  mask_aug     <- NULL
  n_sensors    <- 1L
  cv_eval_mask <- NULL
  dropout_info <- NULL
  n_original   <- nrow(joined_df)

  if (isTRUE(fam$masked)) {
    aug <- spectral_dropout(
      X_proc,
      wavegrid = fam$wavegrid,
      mode     = "profiles",
      profiles = fam$mask_profiles,
      seed     = RUN_SEED
    )
    n_sensors      <- aug$n_copies
    X_proc         <- aug$X
    mask_aug       <- aug$mask
    dataset_labels <- rep(dataset_labels, times = n_sensors)
    sample_ids     <- rep(sample_ids,     times = n_sensors)

    # CV evaluates on full-spectrum original copies only
    cv_eval_mask <- rep(FALSE, nrow(X_proc))
    cv_eval_mask[seq_len(n_original)] <- TRUE

    dropout_info <- list(
      mode = "profiles", n_copies = n_sensors,
      profiles = aug$profiles
    )
    message("  Spectral dropout: ", nrow(aug$profiles),
            " profiles + full-spectrum (", n_sensors, "\u00d7 augmentation)")
    message("  Profiles (nm): ",
            paste(paste0(aug$profiles$wl_lo, "-", aug$profiles$wl_hi),
                  collapse = ", "))
    message("  Augmented corpus: ", nrow(X_proc), " rows")
    message("  CV eval: full-spectrum originals only (n=", n_original, ")")
  }

  # ---- Parallel or sequential training -------------------------------------
  summary_list <- list()
  use_parallel <- PARALLEL && length(props) > 1L

  if (use_parallel) {
    tmp_rds <- tempfile("ossl_worker_", tmpdir = tempdir(), fileext = ".rds")
    saveRDS(list(
      X_proc         = X_proc,
      mask           = mask_aug,
      n_sensors      = n_sensors,
      y_df           = joined_df[, c("Soil_ID", fam$properties), drop = FALSE],
      dataset_labels = dataset_labels,
      sample_ids     = sample_ids,
      cv_eval_mask   = cv_eval_mask,
      family_id      = fam$id,
      family_config  = fam
    ), tmp_rds)
    on.exit(unlink(tmp_rds), add = TRUE)

    # Pre-filter props with enough data
    props_to_run <- Filter(function(p) {
      y    <- suppressWarnings(as.numeric(joined_df[[p]]))
      if (n_sensors > 1L) y <- rep(y, times = n_sensors)
      keep <- is.finite(y) & is.finite(rowSums(X_proc))
      n_ok <- sum(keep)
      if (n_ok < MIN_N) {
        message(sprintf("  [%-20s]  SKIP  (n=%d)", p, n_ok))
        return(FALSE)
      }
      TRUE
    }, props)

    par_results <- dispatch_parallel(
      props           = props_to_run,
      worker_data_rds = tmp_rds,
      family_id       = fam$id,
      max_workers     = N_WORKERS,
      out_dir         = OUT_DIR,
      k               = ACTUAL_K,  epochs    = EPOCHS,
      batch           = BATCH,     pat_es    = PAT_ES,
      pat_lr          = PAT_LR,    latent    = LATENT,
      kl_beta         = KL_BETA,   kl_warmup = KL_WARMUP,
      cv_repeats      = CV_REPEATS,
      seed            = RUN_SEED,
      script_dir      = getwd()
    )

    if (is.null(par_results)) {
      use_parallel <- FALSE
    } else {
      for (prop in names(par_results)) {
        r <- par_results[[prop]]
        if (!is.null(r)) {
          cv <- r$cv
          summary_list[[prop]] <- list(
            property  = prop,
            label     = property_label(prop),
            n         = r$n,
            R2_mean   = cv$R2_mean,   R2_sd   = cv$R2_sd,
            RMSE_mean = cv$RMSE_mean, RMSE_sd = cv$RMSE_sd,
            Bias_mean = cv$Bias_mean, Bias_sd = cv$Bias_sd,
            RPIQ_mean = cv$RPIQ_mean, RPIQ_sd = cv$RPIQ_sd,
            CCC_mean  = cv$CCC_mean,  CCC_sd  = cv$CCC_sd,
            folds     = cv$fold_summary
          )
        }
      }
    }
  }

  # Sequential (single prop, --sequential, or callr unavailable)
  if (!use_parallel) {
    for (prop in props) {
      y    <- suppressWarnings(as.numeric(joined_df[[prop]]))
      if (n_sensors > 1L) y <- rep(y, times = n_sensors)
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
          sample_id   = sample_ids[keep],
          cv_eval_mask = if (!is.null(cv_eval_mask)) cv_eval_mask[keep] else NULL,
          family_id   = fam$id,
          prop        = prop,
          k           = ACTUAL_K,  epochs    = EPOCHS,
          batch       = BATCH,     latent    = LATENT,
          pat_es      = PAT_ES,    pat_lr    = PAT_LR,
          kl_beta     = KL_BETA,   kl_warmup = KL_WARMUP,
          cv_repeats  = CV_REPEATS,
          seed        = RUN_SEED,
          out_dir     = OUT_DIR,
          min_n       = MIN_N,
          mask        = if (!is.null(mask_aug)) mask_aug[keep, , drop = FALSE] else NULL
        ),
        error = function(e) {
          message("    ERROR: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(result)) {
        cv <- result$cv
        message("    ", format_cv_summary(cv))
        for (line in format_cv_fold_lines(cv)) message(line)
        summary_list[[prop]] <- list(
          property  = prop,
          label     = property_label(prop),
          n         = result$n,
          R2_mean   = cv$R2_mean,   R2_sd   = cv$R2_sd,
          RMSE_mean = cv$RMSE_mean, RMSE_sd = cv$RMSE_sd,
          Bias_mean = cv$Bias_mean, Bias_sd = cv$Bias_sd,
          RPIQ_mean = cv$RPIQ_mean, RPIQ_sd = cv$RPIQ_sd,
          CCC_mean  = cv$CCC_mean,  CCC_sd  = cv$CCC_sd,
          folds     = cv$fold_summary
        )
      }

      if (requireNamespace("keras3", quietly = TRUE))
        keras3::clear_session()
      gc(verbose = FALSE)
    }
  }

  # ---- Save metrics summary -------------------------------------------------
  summary_path <- file.path(OUT_DIR, fam$id, "metrics_summary.json")
  dir_create(dirname(summary_path))
  jsonlite::write_json(list(
    family       = fam$id,
    ossl_version = fam$ossl_version,
    pipeline     = pipeline_label,
    masked       = isTRUE(fam$masked),
    dropout      = dropout_info,
    k_folds      = ACTUAL_K,
    cv_repeats   = CV_REPEATS,
    kl_beta      = KL_BETA,
    latent_dim   = LATENT,
    epochs       = EPOCHS,
    seed         = RUN_SEED,
    strategy     = "Grouped stratified K-fold by Soil_ID",
    metrics      = c("RMSE", "Bias", "R2", "RPIQ", "CCC"),
    references   = list(
      RPIQ = "Bellon-Maurel et al. (2010)",
      CCC  = "Lin (1989)",
      OSSL = "Safanelli et al. (2023)"
    ),
    properties   = summary_list
  ), summary_path, auto_unbox = TRUE, pretty = TRUE)

  message("\nSummary saved: ", summary_path)
  message("Done: ", fam$id, " (", length(summary_list), " properties)\n")
}

message(strrep("=", 60))
message("All families complete.")
message(strrep("=", 60))
