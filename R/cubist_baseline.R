# R/cubist_baseline.R - Safanelli 2025 reference model: Cubist + PCA(120) + SNV
#
# Faithful replication of the OSSL modelling pipeline described in
#   Safanelli et al. (2025) PLoS ONE 20(1):e0296545
#
# Pipeline per property:
#   1. SNV preprocessing (already applied upstream by apply_pipeline)
#   2. PCA decomposition keeping the first 120 principal components
#   3. log1p target transform for skewed properties (Safanelli list)
#   4. Cubist fit on the 120 PC scores
#   5. Inner 5-fold CV grid search over committees in {1, 5, 10, 15, 20}
#   6. Save: PCA loadings + means, Cubist model, scaler, Q-statistic UCL,
#            conformal calibration scores, test metrics

# ---- PCA helpers ---------------------------------------------------------

#' Fit a PCA on a spectral matrix and keep the first n_components scores
#'
#' Uses \code{irlba::prcomp_irlba()} when available (randomised truncated SVD —
#' much less memory and faster for rank-120 on large matrices like 130k × 1076)
#' and falls back to \code{stats::prcomp()} otherwise. Returns the fitted
#' object plus centred/scaled metadata for projecting new spectra and computing
#' Q-statistic residuals.
#'
#' @param X Numeric matrix (rows = samples, cols = preprocessed bands)
#' @param n_components Number of PCs to retain (default 120, per Safanelli)
#' @return List with \code{rotation}, \code{center}, \code{scale},
#'   \code{n_components}, \code{scores}, \code{var_explained}
#' @keywords internal
.fit_pca <- function(X, n_components = 120L) {
  X <- as.matrix(X)
  if (ncol(X) < n_components)
    n_components <- ncol(X)

  # Use randomised truncated SVD (irlba) when available: ~10× less memory and
  # faster than full SVD for the first 120 PCs of a 130k × 1076 matrix.
  use_irlba <- requireNamespace("irlba", quietly = TRUE) && nrow(X) > 5000L

  if (use_irlba) {
    # prcomp_irlba returns the same interface as prcomp but uses the IRLBA
    # algorithm; n = number of PCs to compute.
    pc   <- irlba::prcomp_irlba(X, n = n_components, center = TRUE,
                                 scale. = FALSE, retx = TRUE)
    # var_explained: sdev^2 / total variance (estimated from the diagonal of
    # X^T X / (n-1) since we don't compute the full spectrum)
    total_var <- sum(apply(X, 2, var, na.rm = TRUE))
    var_exp   <- (pc$sdev^2) / total_var
  } else {
    pc      <- stats::prcomp(X, center = TRUE, scale. = FALSE,
                             rank. = n_components)
    var_exp <- (pc$sdev[seq_len(n_components)]^2) / sum(pc$sdev^2)
  }

  list(
    rotation      = pc$rotation[, seq_len(n_components), drop = FALSE],
    center        = pc$center,
    scale         = if (is.null(pc$scale) || isFALSE(pc$scale))
                       rep(1, ncol(X)) else pc$scale,
    n_components  = n_components,
    scores        = pc$x[, seq_len(n_components), drop = FALSE],
    var_explained = var_exp[seq_len(n_components)]
  )
}

#' Project new spectra into a fitted PCA space
#'
#' @param pca Object returned by \code{.fit_pca()}
#' @param X_new Numeric matrix in the same band order as the fitting set
#' @return Matrix of PC scores (n_new x n_components)
#' @keywords internal
.project_pca <- function(pca, X_new) {
  X_new <- as.matrix(X_new)
  Xc    <- sweep(X_new, 2, pca$center, "-")
  if (!all(pca$scale == 1)) Xc <- sweep(Xc, 2, pca$scale, "/")
  Xc %*% pca$rotation
}

#' Reconstruct spectra from the first n PC scores (for Q-statistic)
#' @keywords internal
.reconstruct_pca <- function(pca, scores) {
  X_hat <- scores %*% t(pca$rotation)
  X_hat <- sweep(X_hat, 2, pca$center, "+")
  X_hat
}

# ---- Cubist trainer with inner-CV committee tuning -----------------------

# Default Cubist committee grid per Safanelli (Methods, p.9)
.CUBIST_COMMITTEE_GRID <- c(1L, 5L, 10L, 15L, 20L)

#' Inner 5-fold CV to choose the best Cubist committees value
#'
#' Per Safanelli: hyperparameter search on a SUBSET of ~2000 samples for
#' speed. We mirror that with \code{inner_subset_n = 2000}.
#'
#' @param T_scores Matrix of PC scores (rows = samples, cols = PCs)
#' @param y Numeric target vector
#' @param committees_grid Integer vector of committee values to try
#' @param k_inner Inner-CV folds (default 5, per Safanelli)
#' @param inner_subset_n Subsample size for tuning (default 2000)
#' @param seed RNG seed
#' @return List with \code{best_committees} and per-config RMSE means
#' @keywords internal
.tune_cubist_committees <- function(T_scores, y,
                                    committees_grid = .CUBIST_COMMITTEE_GRID,
                                    k_inner         = 5L,
                                    inner_subset_n  = 2000L,
                                    seed            = 42L) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required for the Cubist baseline. ",
         "Install with: install.packages('Cubist')")

  set.seed(seed)
  n <- nrow(T_scores)
  if (n > inner_subset_n) {
    sub_idx <- sample.int(n, inner_subset_n)
    T_in    <- T_scores[sub_idx, , drop = FALSE]
    y_in    <- y[sub_idx]
  } else {
    T_in <- T_scores; y_in <- y
  }

  fold_id <- sample(rep(seq_len(k_inner), length.out = nrow(T_in)))
  rmse_mat <- matrix(NA_real_,
                     nrow = length(committees_grid), ncol = k_inner,
                     dimnames = list(as.character(committees_grid), NULL))

  for (k in seq_len(k_inner)) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    for (ci in seq_along(committees_grid)) {
      cm <- committees_grid[ci]
      mdl <- Cubist::cubist(x = T_in[tr, , drop = FALSE],
                            y = y_in[tr],
                            committees = cm)
      yh  <- predict(mdl, T_in[te, , drop = FALSE])
      rmse_mat[ci, k] <- sqrt(mean((y_in[te] - yh)^2, na.rm = TRUE))
    }
  }
  mean_rmse <- rowMeans(rmse_mat, na.rm = TRUE)
  best <- committees_grid[which.min(mean_rmse)]
  list(best_committees = as.integer(best),
       mean_rmse       = setNames(as.numeric(mean_rmse),
                                  as.character(committees_grid)))
}

# ---- Single-property Cubist+PCA trainer ----------------------------------

#' Train one Cubist + PCA(120) model for a single soil property
#'
#' Faithful to Safanelli 2025: PCA(120) on SNV-preprocessed spectra,
#' optional log1p target transform, Cubist with inner-CV committee tuning,
#' Q-statistic UCL on training residuals, normalised-conformal calibration
#' for prediction intervals.
#'
#' @param X Numeric matrix of preprocessed (SNV) spectra
#' @param y Numeric target vector
#' @param family_id Family ID (e.g. "OSSL_VisNIR_Cubist")
#' @param prop Soil property name
#' @param out_dir Output root (default "models")
#' @param n_components PCA components (default 120)
#' @param committees_grid Cubist committee values to tune
#' @param k_inner Inner-CV folds for tuning (default 5)
#' @param inner_subset_n Subsample for tuning (default 2000)
#' @param min_n Minimum samples to attempt training (default 50)
#' @return Invisible TRUE on success, NULL on skip
#' @export
train_cubist_pca <- function(X, y, family_id, prop,
                             out_dir          = "models",
                             n_components     = 120L,
                             committees_grid  = .CUBIST_COMMITTEE_GRID,
                             k_inner          = 5L,
                             inner_subset_n   = 2000L,
                             min_n            = 50L) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required. install.packages('Cubist')")

  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  # 80/20 train/test split (matches Safanelli's 10-fold structure: each fold
  # uses 90/10, so 80/20 is the closest single-split equivalent that keeps
  # the test set sized for sane metric estimation).
  set.seed(42L)
  idx  <- sample.int(n)
  n_tr <- floor(0.80 * n)
  X_tr <- X[idx[seq_len(n_tr)],          , drop = FALSE]; y_tr <- y[idx[seq_len(n_tr)]]
  X_te <- X[idx[-seq_len(n_tr)],         , drop = FALSE]; y_te <- y[idx[-seq_len(n_tr)]]

  # 1. PCA on training set
  pca <- .fit_pca(X_tr, n_components = n_components)

  # 2. Inner-CV committee tuning
  tune <- .tune_cubist_committees(pca$scores, y_tr,
                                  committees_grid = committees_grid,
                                  k_inner         = k_inner,
                                  inner_subset_n  = inner_subset_n)

  # 3. Final fit on full training set with best committees
  mdl <- Cubist::cubist(x = pca$scores, y = y_tr,
                        committees = tune$best_committees)

  # 4. Test set evaluation — metrics in ORIGINAL scale (back-transform log1p if
  #    needed so they are comparable to Safanelli 2025 Tables 4/5 which report
  #    original-unit RMSE, R², CCC, RPIQ for all properties).
  T_te    <- .project_pca(pca, X_te)
  yhat_te <- predict(mdl, T_te)
  if (prop %in% ossl_l1_log_transform) {
    mets <- metrics_from_y(expm1(y_te), expm1(yhat_te))
  } else {
    mets <- metrics_from_y(y_te, yhat_te)
  }

  # 5. Q-statistic UCL from training residuals
  X_tr_hat <- .reconstruct_pca(pca, pca$scores)
  E_tr     <- X_tr - X_tr_hat
  ucl_q    <- q_statistic_ucl(E_tr, alpha = 0.01)   # 99% conf, Safanelli

  # 6. Normalised conformal — Safanelli-faithful: residuals come from
  #    INTERNAL 5-fold CV on the training set, NOT from a held-out calib set.
  conf <- train_normalized_conformal(
    pca            = pca,
    response_model = mdl,
    X_train        = X_tr,
    y_train        = y_tr,
    confidence     = 0.68,
    committees     = tune$best_committees,
    k_cv           = 5L
  )

  # Calibration-quality residuals (used for q90/q95 storage only)
  yhat_tr <- as.numeric(predict(mdl, pca$scores))
  res_ca  <- abs(y_tr - yhat_tr)

  # 8. Save artefacts
  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  saveRDS(list(model           = mdl,
               pca             = pca,
               best_committees = tune$best_committees),
          file = paste0(base, "_cubist.rds"))

  # Error model serialised separately (Cubist object isn't JSON-serialisable)
  saveRDS(conf$error_model, file = paste0(base, "_cubist_error.rds"))

  # 9. Metrics + applicability + conformal info (scalars only in JSON)
  fr <- feature_minmax(X)
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = as.numeric(quantile(res_ca, 0.90, na.rm = TRUE)),
                                      q95 = as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))),
                    feat_range = fr,
                    latent     = list(
                      method        = "q_statistic",
                      ucl_q         = ucl_q,
                      n_components  = n_components,
                      var_explained = sum(pca$var_explained),
                      committees    = tune$best_committees,
                      conformal     = list(
                        alpha_q    = conf$alpha_q,
                        confidence = conf$confidence
                      )
                    ))

  message(sprintf("  [%s/%s] RMSE=%.3f Bias=%.3f R2=%.3f RPIQ=%.2f CCC=%.3f committees=%d n=%d",
                  family_id, prop,
                  mets$RMSE, mets$Bias, mets$R2, mets$RPIQ, mets$CCC,
                  tune$best_committees, n))
  invisible(TRUE)
}

# ---- Cubist+PCA loader ---------------------------------------------------

#' Load a saved Cubist + PCA bundle for a property
#'
#' @param family_id Family ID
#' @param prop Soil property name
#' @param model_dir Root model directory (default "models")
#' @return List with \code{model} (Cubist), \code{pca}, \code{best_committees}
#' @export
load_cubist_pca <- function(family_id, prop, model_dir = "models") {
  base <- file.path(model_dir, family_id, "models", prop)
  rds  <- paste0(base, "_cubist.rds")
  if (!file.exists(rds))
    stop("Cubist bundle not found: ", rds)
  readRDS(rds)
}

#' Predict with a Cubist + PCA bundle
#'
#' @param bundle List from \code{load_cubist_pca()}
#' @param X_new Preprocessed spectral matrix in the same band order as training
#' @return Numeric vector of predictions
#' @export
predict_cubist_pca <- function(bundle, X_new) {
  T_new <- .project_pca(bundle$pca, X_new)
  as.numeric(predict(bundle$model, T_new))
}
