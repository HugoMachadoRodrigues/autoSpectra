# R/train_cv.R - 10-fold cross-validation with refitting (Safanelli 2025)
#
# Per Safanelli (Methods, p.10): "The evaluation of the models was performed
# with a 10-fold external cross-validation (outer) of the tuned models using
# RMSE, mean error (bias), R-squared (R2), Lin's CCC, and RPIQ."
#
# Strategy:
#   - 10 random folds across the full OSSL training corpus
#   - For each fold: hold out test set, refit response model on remaining 9
#   - log1p target transform applied for skewed properties
#   - Return per-fold metrics (mean +/- SD reported in metrics_summary.json)
#
# An optional --lodo flag (Leave-One-Dataset-Out) is also supported for
# stricter cross-instrument evaluation, where each fold withholds an entire
# OSSL contributing dataset. This is a more honest cross-sensor estimate
# than random k-fold.

# ---- Fold construction ---------------------------------------------------

#' Build random 10-fold indices (Safanelli default)
#'
#' @param n Number of samples
#' @param k Number of folds (default 10)
#' @param seed RNG seed
#' @return List of length k; each element is integer vector of test indices
#' @export
random_kfold <- function(n, k = 10L, seed = 42L) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(k), length.out = n))
  lapply(seq_len(k), function(fold) which(fold_id == fold))
}

#' Build Leave-One-Dataset-Out (LODO) k-fold indices
#'
#' Each fold withholds all samples from one OSSL contributing dataset.
#' Datasets are randomly assigned to folds when there are more datasets
#' than folds.
#'
#' @param dataset_col Character vector of dataset labels (one per sample)
#' @param k Number of folds (default 10)
#' @param seed RNG seed
#' @return List of length k of integer test-index vectors
#' @export
stratified_kfold <- function(dataset_col, k = 10L, seed = 42L) {
  set.seed(seed)
  datasets <- unique(dataset_col)
  ds_order <- sample(datasets)
  ds_fold  <- setNames(((seq_along(ds_order) - 1L) %% k) + 1L, ds_order)
  folds    <- lapply(seq_len(k),
                     function(f) which(ds_fold[dataset_col] == f))
  folds
}

# ---- Per-property k-fold CV (Safanelli outer loop) -----------------------

#' Perform 10-fold CV with refitting for one (family, property) combination
#'
#' Dispatches to either VAE or Cubist+PCA training on each fold based on the
#' family's algorithm.
#'
#' @param X Numeric matrix of preprocessed (SNV) spectra
#' @param y Numeric target vector
#' @param family_id Family identifier
#' @param prop Soil property name
#' @param folds List of integer test-index vectors (from \code{random_kfold()}
#'   or \code{stratified_kfold()})
#' @param fam Family list from \code{model_registry}
#' @param ... Extra args forwarded to per-fold trainers
#' @return List with per-fold metrics and aggregated mean/SD
#' @export
cv_train_property <- function(X, y, family_id, prop, folds, fam, ...) {
  algo  <- fam$algorithm %||% "vae"
  k     <- length(folds)
  per_fold <- vector("list", k)

  for (i in seq_len(k)) {
    te <- folds[[i]]
    tr <- setdiff(seq_along(y), te)
    if (length(tr) < 50L || length(te) < 5L) {
      per_fold[[i]] <- NULL
      next
    }

    if (identical(algo, "vae_masked")) {
      mets <- .cv_one_fold_vae_masked(X[tr, ], y[tr], X[te, ], y[te],
                                      family_id, prop, fam, ...)
    } else if (identical(algo, "vae")) {
      mets <- .cv_one_fold_vae(X[tr, ], y[tr], X[te, ], y[te],
                               family_id, prop, ...)
    } else if (identical(algo, "cubist_pca")) {
      mets <- .cv_one_fold_cubist(X[tr, ], y[tr], X[te, ], y[te],
                                  prop = prop, ...)
    } else {
      stop("Unknown algorithm: ", algo)
    }
    per_fold[[i]] <- mets
  }

  ok        <- !vapply(per_fold, is.null, logical(1))
  fold_mets <- per_fold[ok]
  if (length(fold_mets) == 0L) return(NULL)

  agg <- .aggregate_fold_metrics(fold_mets)
  list(per_fold = fold_mets, mean_sd = agg)
}

# ---- Per-fold workers ----------------------------------------------------

.cv_one_fold_cubist <- function(X_tr, y_tr, X_te, y_te,
                                prop            = NULL,
                                n_components    = 120L,
                                committees_grid = c(1L, 5L, 10L, 15L, 20L),
                                k_inner         = 5L,
                                inner_subset_n  = 2000L,
                                ...) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required.")
  pca   <- .fit_pca(X_tr, n_components = n_components)
  tune  <- .tune_cubist_committees(pca$scores, y_tr,
                                   committees_grid = committees_grid,
                                   k_inner         = k_inner,
                                   inner_subset_n  = inner_subset_n)
  mdl   <- Cubist::cubist(x = pca$scores, y = y_tr,
                          committees = tune$best_committees)
  T_te  <- .project_pca(pca, X_te)
  yhat  <- as.numeric(predict(mdl, T_te))
  # Back-transform log1p for original-scale metrics (Safanelli comparability)
  mets <- metrics_original_scale(y_te, yhat, prop, bounds = target_bounds(y_tr))
  mets$committees <- tune$best_committees
  mets$n_test     <- length(y_te)
  mets
}

.cv_one_fold_vae <- function(X_tr, y_tr, X_te, y_te,
                             family_id, prop,
                             epochs           = 80L,
                             batch_size       = 32L,
                             val_split        = 0.15,
                             patience_es      = 10L,
                             patience_rl      = 5L,
                             latent_dim       = 16L,
                             kl_beta          = 1e-3,
                             kl_warmup_epochs = 30L,
                             ...) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  mean_y <- mean(y_tr, na.rm = TRUE)
  sd_y   <- sd(y_tr, na.rm = TRUE); if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z    <- (y_tr - mean_y) / sd_y

  mdl <- build_soilVAE(ncol(X_tr), latent_dim = latent_dim)
  cbs <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = patience_es, restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = patience_rl, factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup_epochs))
  )
  mdl |> keras3::fit(
    x                = as.matrix(X_tr),
    y                = list(as.matrix(X_tr), y_z),
    epochs           = epochs,
    batch_size       = batch_size,
    validation_split = val_split,
    callbacks        = cbs,
    verbose          = 0
  )
  yhat <- .extract_prediction(mdl, X_te) * sd_y + mean_y
  mets <- metrics_original_scale(y_te, yhat, prop, bounds = target_bounds(y_tr))
  mets$n_test <- length(y_te)
  mets
}

.cv_one_fold_vae_masked <- function(X_tr, y_tr, X_te, y_te,
                                    family_id, prop, fam,
                                    epochs           = 80L,
                                    batch_size       = 64L,
                                    val_split        = 0.15,
                                    patience_es      = 10L,
                                    patience_rl      = 5L,
                                    latent_dim       = 16L,
                                    kl_beta          = 1e-3,
                                    kl_warmup_epochs = 30L,
                                    n_aug_copies     = 2L,
                                    ...) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  mean_y <- mean(y_tr, na.rm = TRUE)
  sd_y   <- sd(y_tr, na.rm = TRUE); if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z    <- (y_tr - mean_y) / sd_y

  aug   <- spectral_dropout(X_tr, wavegrid = fam$wavegrid,
                            sensor_type = fam$sensor_type,
                            n_copies = n_aug_copies, seed = 42L)
  y_aug <- rep(y_z, times = n_aug_copies + 1L)

  X_in  <- .pack_masked_input(aug$X, aug$mask)
  Y_rec <- .pack_masked_target(aug$X, aug$mask)

  d_spec <- ncol(X_tr)
  mdl <- build_maskedVAE(d_spec, latent_dim = latent_dim)
  cbs <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = patience_es, restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = patience_rl, factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup_epochs))
  )
  mdl |> keras3::fit(
    x                = X_in,
    y                = list(Y_rec, y_aug),
    epochs           = epochs,
    batch_size       = batch_size,
    validation_split = val_split,
    callbacks        = cbs,
    verbose          = 0
  )
  X_te_in <- .pack_masked_input(X_te, matrix(1L, nrow(X_te), d_spec))
  yhat    <- .extract_prediction(mdl, X_te_in) * sd_y + mean_y
  mets <- metrics_original_scale(y_te, yhat, prop, bounds = target_bounds(y_tr))
  mets$n_test <- length(y_te)
  mets
}

# ---- Aggregation ---------------------------------------------------------

.aggregate_fold_metrics <- function(fold_mets) {
  keys <- c("RMSE", "Bias", "R2", "RPIQ", "CCC")
  out  <- list()
  for (k in keys) {
    vals <- vapply(fold_mets,
                   function(m) suppressWarnings(as.numeric(m[[k]])),
                   numeric(1))
    out[[paste0(k, "_mean")]] <- mean(vals, na.rm = TRUE)
    out[[paste0(k, "_sd")]]   <- stats::sd(vals, na.rm = TRUE)
  }
  out$n_folds <- length(fold_mets)
  out$n_test  <- sum(vapply(fold_mets,
                            function(m) m$n_test %||% NA_integer_,
                            integer(1)), na.rm = TRUE)
  out
}
