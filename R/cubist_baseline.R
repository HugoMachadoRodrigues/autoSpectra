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
#' @param scale Divide each band by its standard deviation before the PCA, as
#'   the OSSL models do (\code{prcomp(center = TRUE, scale = TRUE)}); bands
#'   with no spread keep a divisor of 1
#' @return List with \code{rotation}, \code{center}, \code{scale},
#'   \code{n_components}, \code{scores}, \code{var_explained}
#' @keywords internal
.fit_pca <- function(X, n_components = 120L, scale = FALSE) {
  X <- as.matrix(X)
  if (ncol(X) < n_components)
    n_components <- ncol(X)
  col_var <- apply(X, 2, stats::var, na.rm = TRUE)
  sds <- FALSE
  if (isTRUE(scale)) {
    sds <- sqrt(col_var)
    sds[!is.finite(sds) | sds <= 0] <- 1
  }

  # Use randomised truncated SVD (irlba) when available: ~10× less memory and
  # faster than full SVD for the first 120 PCs of a 130k × 1076 matrix.
  use_irlba <- requireNamespace("irlba", quietly = TRUE) && nrow(X) > 5000L

  if (use_irlba) {
    # prcomp_irlba returns the same interface as prcomp but uses the IRLBA
    # algorithm; n = number of PCs to compute.
    pc   <- irlba::prcomp_irlba(X, n = n_components, center = TRUE,
                                 scale. = sds, retx = TRUE)
    # var_explained: sdev^2 / total variance (the trace of the covariance of
    # the centred, possibly scaled, matrix; the full spectrum is not computed)
    total_var <- sum(if (isFALSE(sds)) col_var else col_var / sds^2)
    var_exp   <- (pc$sdev^2) / total_var
  } else {
    pc      <- stats::prcomp(X, center = TRUE, scale. = sds,
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

# ---- Redistribution guard ------------------------------------------------

#' Strip embedded training data from a fitted Cubist model
#'
#' \code{Cubist::cubist()} keeps a verbatim copy of the whole training matrix
#' (all predictor columns plus the response) in \code{$data}, as a CSV text
#' blob. Because we train on PCA scores of OSSL spectra, that blob would let
#' anyone reconstruct the OSSL spectra (via \code{.reconstruct_pca()}) and read
#' the OSSL lab values straight out of a published model file — i.e. it would
#' turn our model archive into an unlicensed redistribution of the OSSL
#' database. It is also pure dead weight: predictions are bit-identical
#' without it.
#'
#' Always call this before serialising a Cubist model.
#'
#' @param mdl A model returned by \code{Cubist::cubist()}
#' @return The same model with \code{$data} emptied
#' @keywords internal
.strip_cubist_data <- function(mdl) {
  if (is.null(mdl)) return(mdl)
  mdl$data <- ""
  mdl
}

#' Drop the per-sample PC scores from a fitted PCA before serialisation
#'
#' \code{pca$scores} holds the 120 PC scores of every OSSL training sample and
#' \code{.reconstruct_pca()} inverts them back to the input spectra at ~99.98%
#' of the variance. Only the projection basis (rotation/center/scale) is needed
#' at inference time — \code{.project_pca()} never reads \code{scores}.
#'
#' @param pca A list from \code{.fit_pca()}
#' @return The same list with \code{scores} replaced by a row count
#' @keywords internal
.strip_pca_scores <- function(pca) {
  n_train      <- if (is.null(pca$scores)) NA_integer_ else nrow(pca$scores)
  pca$scores   <- NULL
  pca$n_train  <- n_train    # provenance only
  pca
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
#' @param groups Optional group label per row of \code{T_scores}; rows sharing
#'   one are kept in the same inner fold
#' @return List with \code{best_committees} and per-config RMSE means
#' @keywords internal
.tune_cubist_committees <- function(T_scores, y,
                                    committees_grid = .CUBIST_COMMITTEE_GRID,
                                    k_inner         = 5L,
                                    inner_subset_n  = 2000L,
                                    seed            = 42L,
                                    groups          = NULL) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required for the Cubist baseline. ",
         "Install with: install.packages('Cubist')")

  set.seed(seed)
  n <- nrow(T_scores)
  if (is.null(groups)) groups <- seq_len(n)
  sub_idx <- if (n > inner_subset_n) sample.int(n, inner_subset_n) else seq_len(n)
  T_in <- T_scores[sub_idx, , drop = FALSE]
  y_in <- y[sub_idx]

  folds   <- random_kfold(nrow(T_in), k = k_inner, seed = seed,
                          groups = groups[sub_idx])
  fold_id <- integer(nrow(T_in))
  for (k in seq_along(folds)) fold_id[folds[[k]]] <- k
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
#' @param groups Optional group label per row; rows sharing one stay on the
#'   same side of the train/test split and of every internal CV fold
#' @return Invisible TRUE on success, NULL on skip
#' @export
train_cubist_pca <- function(X, y, family_id, prop,
                             out_dir          = "models",
                             n_components     = 120L,
                             committees_grid  = .CUBIST_COMMITTEE_GRID,
                             k_inner          = 5L,
                             inner_subset_n   = 2000L,
                             min_n            = 50L,
                             groups           = NULL) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required. install.packages('Cubist')")

  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  # 80/20 train/test split (matches Safanelli's 10-fold structure: each fold
  # uses 90/10, so 80/20 is the closest single-split equivalent that keeps
  # the test set sized for sane metric estimation). Rows sharing a group (repeat
  # visits of one LUCAS point) stay on the same side.
  if (is.null(groups)) groups <- seq_len(n)
  sp   <- .grouped_partition(n, c(0.80, 0.20), groups, seed = 42L)
  X_tr <- X[sp[[1]], , drop = FALSE]; y_tr <- y[sp[[1]]]; g_tr <- groups[sp[[1]]]
  X_te <- X[sp[[2]], , drop = FALSE]; y_te <- y[sp[[2]]]
  # The held-out 20% is halved for the prediction intervals: one half
  # calibrates them, the other measures their coverage. Point metrics use all
  # of it, since neither half is seen by the response model.
  half   <- .grouped_partition(length(y_te), c(0.5, 0.5), groups[sp[[2]]], seed = 44L)
  cal_ix <- half[[1]]; eval_ix <- half[[2]]

  # 1. PCA on training set
  pca <- .fit_pca(X_tr, n_components = n_components)

  # 2. Inner-CV committee tuning
  tune <- .tune_cubist_committees(pca$scores, y_tr,
                                  committees_grid = committees_grid,
                                  k_inner         = k_inner,
                                  inner_subset_n  = inner_subset_n,
                                  groups          = g_tr)

  # 3. Final fit on full training set with best committees
  mdl <- Cubist::cubist(x = pca$scores, y = y_tr,
                        committees = tune$best_committees)

  # 4. Test set evaluation — metrics in ORIGINAL units, so they are comparable
  #    to Safanelli 2025 Tables 4/5. The back-transform is bounded by the
  #    training range: see R/target_transform.R for why an unbounded expm1()
  #    destroys log1p properties under domain shift.
  bounds  <- target_bounds(y_tr)
  T_te    <- .project_pca(pca, X_te)
  yhat_te <- predict(mdl, T_te)
  mets    <- metrics_original_scale(y_te, yhat_te, prop, bounds = bounds)

  # 5. Applicability-domain limit. Safanelli uses the Jackson-Mudholkar UCL,
  #    but on OSSL it is badly calibrated: residuals are far from the normal
  #    approximation it assumes, and 9-37% of TRAINING spectra exceeded the
  #    nominal 1% limit. The limit used is therefore the empirical 99th
  #    percentile of the training Q; the J-M value is kept for reference.
  E_tr     <- X_tr - .reconstruct_pca(pca, pca$scores)
  Q_tr     <- q_statistic(E_tr)
  ucl_q    <- as.numeric(stats::quantile(Q_tr, 0.99, names = FALSE))
  ucl_q_jm <- q_statistic_ucl(E_tr, alpha = 0.01)
  rm(E_tr)

  # 6. Normalised conformal — residuals from internal CV on the training set
  #    (Safanelli), with a separate half for the non-conformity scores.
  conf <- train_normalized_conformal(
    pca            = pca,
    response_model = mdl,
    X_train        = X_tr,
    y_train        = y_tr,
    confidence     = 0.68,
    committees     = tune$best_committees,
    k_cv           = 5L,
    groups         = g_tr,
    X_cal          = X_te[cal_ix, , drop = FALSE],
    y_cal          = y_te[cal_ix]
  )

  # 7. Held-out diagnostics on the evaluation half: how often spectra are
  #    flagged out of domain, and how often the intervals cover the observed
  #    value. Both should be near nominal.
  Q_te      <- predict_q_statistic(pca, X_te[eval_ix, , drop = FALSE], ucl_q)
  delta_ev  <- pmax(as.numeric(predict(conf$error_model, T_te[eval_ix, , drop = FALSE])), 0)
  abs_err   <- abs(y_te[eval_ix] - yhat_te[eval_ix])
  coverage  <- vapply(c(0.68, 0.90, 0.95), function(lv)
    mean(abs_err <= delta_ev * .conformal_alpha(conf, lv), na.rm = TRUE), numeric(1))

  # 8. Save artefacts
  #
  # REDISTRIBUTION GUARD — both the PCA scores and Cubist's embedded `$data`
  # blob would let anyone reconstruct the OSSL spectra and read the OSSL lab
  # values out of a published model file. Strip both before serialising; see
  # .strip_pca_scores() / .strip_cubist_data(). Neither affects predictions.
  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  saveRDS(list(model           = .strip_cubist_data(mdl),
               pca             = .strip_pca_scores(pca),
               best_committees = tune$best_committees),
          file = paste0(base, "_cubist.rds"))

  # Error model serialised separately (Cubist object isn't JSON-serialisable)
  saveRDS(.strip_cubist_data(conf$error_model),
          file = paste0(base, "_cubist_error.rds"))

  # Cubist needs no y-scaling, but prediction still needs the calibration range
  # to bound the back-transform. Write it to the same _scaler path the VAE
  # families use, so predict_soil() finds it the same way for every family.
  save_scaler(0, 1, path_base = base, bounds = bounds)

  # 9. Metrics + applicability + conformal info (aggregate values only)
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(),
                    latent     = list(
                      method          = "q_statistic",
                      ucl_q           = ucl_q,
                      ucl_q_method    = "empirical_q99_train",
                      ucl_q_jm        = ucl_q_jm,
                      train_frac_above_jm = mean(Q_tr > ucl_q_jm),
                      test_frac_out_of_domain = mean(!Q_te$in_domain),
                      n_components    = n_components,
                      var_explained   = sum(pca$var_explained),
                      committees      = tune$best_committees,
                      n_train         = length(y_tr),
                      n_test          = length(y_te),
                      conformal       = list(
                        alpha_q       = conf$alpha_q,
                        confidence    = conf$confidence,
                        levels        = conf$levels,
                        alpha_grid    = conf$alpha_grid,
                        n_calibration = conf$n_calibration,
                        test_coverage_68 = coverage[1],
                        test_coverage_90 = coverage[2],
                        test_coverage_95 = coverage[3],
                        n_coverage_eval  = length(eval_ix)
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
  # Loading the namespace registers Cubist's predict() method
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required. install.packages('Cubist')")
  T_new <- .project_pca(bundle$pca, X_new)
  as.numeric(predict(bundle$model, T_new))
}
