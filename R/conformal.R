# R/conformal.R - Normalised conformal prediction (Safanelli 2025, eq. 9-10)
#
# Per Safanelli, uncertainty is calibrated by training a SECOND model that
# predicts the absolute error of the response model:
#
#   delta_i = |y_i - y_hat_i|         (cross-validated residual)
#   delta_hat = f(T120)                (error model, same structure as response)
#   alpha_i = |y_i - y_hat_i| / delta_hat_i        (eq. 9: non-conformity score)
#   alpha_q  = quantile(alpha, q)                  (q = confidence, e.g. 0.68)
#
# Prediction interval for a new sample:
#   PI = y_hat +/- delta_hat * alpha_q              (eq. 10)
#
# 68% confidence approximates 1 SD - this is what Safanelli reports.

#' Train a normalised-conformal calibration for a Cubist + PCA model
#'
#' Mirrors the Safanelli 2025 procedure (eqs. 9-10):
#' \enumerate{
#'   \item Cross-validated predictions \eqn{\hat{y}_i^{CV}} on the training
#'         set via internal k-fold CV (Safanelli: "Cross-validated predictions
#'         were also used to estimate the unbiased absolute error", p.10)
#'   \item Compute absolute residuals \eqn{\delta_i = |y_i - \hat{y}_i^{CV}|}
#'   \item Train an error model with the SAME architecture (Cubist on the same
#'         PC scores) that maps PC scores to \eqn{\delta_i}
#'   \item Compute non-conformity scores \eqn{\alpha_i = \delta_i / \hat{\delta}_i}
#'   \item Store the \code{confidence}-th quantile of \eqn{\alpha_i}
#' }
#'
#' @param pca Fitted PCA list (from \code{.fit_pca()})
#' @param response_model Cubist model trained on \code{pca$scores}
#' @param X_train Preprocessed TRAINING spectra used to fit \code{response_model}
#' @param y_train Numeric vector of training targets
#' @param confidence Confidence level (default 0.68 = ~1 SD per Safanelli)
#' @param committees Cubist committees for the error model (default = same as
#'   response model when retrievable, else 5)
#' @param k_cv Folds used to obtain CV residuals on the training set (default 5)
#' @param conformal_subset_n Cap on how many training samples feed the inner CV
#'   that produces the residuals (default 10000). Calibration only needs a stable
#'   quantile, so a subsample avoids refitting the full training set k_cv times.
#' @param seed RNG seed
#' @return List with \code{error_model}, \code{alpha_q}, \code{confidence}
#' @export
train_normalized_conformal <- function(pca, response_model,
                                       X_train, y_train,
                                       confidence    = 0.68,
                                       committees    = NULL,
                                       k_cv          = 5L,
                                       conformal_subset_n = 10000L,
                                       seed          = 42L) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required for normalized conformal calibration.")

  cm <- if (!is.null(committees)) committees
        else if (!is.null(response_model$committees)) response_model$committees
        else 5L

  # 1. Internal k-fold CV to obtain unbiased residuals on training set
  #    (Safanelli p.10: cross-validated predictions are the CONFORMAL source,
  #    not a held-out calibration set).
  #
  #    For very large n, we subsample to keep the inner-CV cost reasonable.
  #    Safanelli reuses outer-loop CV predictions (free); single-split mode
  #    needs a fresh CV, hence the subsample.  10 000 samples gives stable
  #    quantile estimates without the 5x full-train cost.
  set.seed(seed)
  T_full <- .project_pca(pca, X_train)
  n_full <- nrow(T_full)

  if (n_full > conformal_subset_n) {
    sub_idx <- sample.int(n_full, conformal_subset_n)
    T_in    <- T_full[sub_idx, , drop = FALSE]
    y_in    <- y_train[sub_idx]
  } else {
    T_in <- T_full; y_in <- y_train
  }
  n_in    <- nrow(T_in)
  fold_id <- sample(rep(seq_len(k_cv), length.out = n_in))
  yhat_cv <- numeric(n_in)
  for (f in seq_len(k_cv)) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    mdl_f <- Cubist::cubist(x = T_in[tr, , drop = FALSE],
                            y = y_in[tr],
                            committees = cm)
    yhat_cv[te] <- as.numeric(predict(mdl_f, T_in[te, , drop = FALSE]))
  }
  delta <- abs(y_in - yhat_cv)
  delta[!is.finite(delta)] <- NA_real_

  # 2. Error model fit on the SAME subsample's PC scores against CV residuals
  err_mdl <- Cubist::cubist(x = T_in, y = delta, committees = cm)

  # 3. Non-conformity scores on the subsample (uses CV residuals)
  delta_hat <- as.numeric(predict(err_mdl, T_in))
  delta_hat <- pmax(delta_hat, .Machine$double.eps)
  alpha     <- delta / delta_hat
  alpha     <- alpha[is.finite(alpha)]
  alpha_q   <- as.numeric(stats::quantile(alpha, confidence, na.rm = TRUE))

  list(
    error_model        = err_mdl,
    alpha_q            = alpha_q,
    confidence         = confidence,
    cv_folds           = k_cv,
    conformal_subset_n = min(n_full, conformal_subset_n)
  )
}

#' Predict with prediction intervals using a normalised-conformal calibration
#'
#' @param pca Fitted PCA
#' @param response_model Cubist response model
#' @param error_model Cubist error model from \code{train_normalized_conformal()}
#' @param alpha_q Calibrated non-conformity quantile (alpha_q)
#' @param X_new Preprocessed spectra (matrix)
#' @return Data frame with columns: y_pred, delta_hat, pi_lo, pi_hi
#' @export
predict_with_pi <- function(pca, response_model, error_model, alpha_q, X_new) {
  T_new      <- .project_pca(pca, X_new)
  y_pred     <- as.numeric(predict(response_model, T_new))
  delta_hat  <- as.numeric(predict(error_model, T_new))
  delta_hat  <- pmax(delta_hat, 0)
  data.frame(
    y_pred    = y_pred,
    delta_hat = delta_hat,
    pi_lo     = y_pred - delta_hat * alpha_q,
    pi_hi     = y_pred + delta_hat * alpha_q,
    stringsAsFactors = FALSE
  )
}

#' Simple absolute-residual conformal calibration (used by the VAE family)
#'
#' Less sophisticated than the normalised variant, but does not require a
#' second model. Returns a single half-width per confidence level.
#'
#' @param y_true Calibration observed values
#' @param y_pred Calibration predicted values
#' @param confidence Numeric vector of confidence levels (default c(0.68, 0.90, 0.95))
#' @return Named numeric vector of half-widths (one per confidence level)
#' @export
absolute_conformal <- function(y_true, y_pred,
                               confidence = c(0.68, 0.90, 0.95)) {
  res <- abs(y_true - y_pred)
  res <- res[is.finite(res)]
  setNames(as.numeric(stats::quantile(res, confidence, na.rm = TRUE)),
           paste0("q", round(confidence * 100)))
}
