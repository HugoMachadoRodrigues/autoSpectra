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

# Confidence levels whose non-conformity quantile is stored with each model, so
# predict_soil(confidence = ...) can serve any of them. These are aggregate
# statistics of the calibration scores, not per-sample values.
.CONFORMAL_LEVELS <- seq(0.50, 0.99, by = 0.01)

#' Train a normalised-conformal calibration for a Cubist + PCA model
#'
#' Follows the Safanelli 2025 procedure (eqs. 9-10), with the calibration split
#' that makes it a proper split-conformal method:
#' \enumerate{
#'   \item Cross-validated predictions \eqn{\hat{y}_i^{CV}} on (a subsample of)
#'         the training set via internal k-fold CV
#'   \item Absolute residuals \eqn{\delta_i = |y_i - \hat{y}_i^{CV}|}
#'   \item The subsample is halved. An error model (Cubist on the same PC
#'         scores) learns \eqn{\delta} on one half...
#'   \item ...and the non-conformity scores
#'         \eqn{\alpha_i = \delta_i / \hat{\delta}_i} are computed on the other.
#'         Scoring the rows the error model was fitted on understates
#'         \eqn{\hat{\delta}} errors and gives intervals that over-cover.
#'   \item The finite-sample-corrected quantile of \eqn{\alpha} is stored for
#'         every level in \code{levels}
#' }
#'
#' @param pca Fitted PCA list (from \code{.fit_pca()})
#' @param response_model Cubist model trained on \code{pca$scores}
#' @param X_train Preprocessed TRAINING spectra used to fit \code{response_model}
#' @param y_train Numeric vector of training targets
#' @param confidence Default confidence level (0.68 = ~1 SD per Safanelli)
#' @param committees Cubist committees for the error model (default = same as
#'   response model when retrievable, else 5)
#' @param k_cv Folds used to obtain CV residuals on the training set (default 5)
#' @param conformal_subset_n Cap on how many training samples feed the inner CV
#'   that produces the residuals (default 10000). Calibration only needs a stable
#'   quantile, so a subsample avoids refitting the full training set k_cv times.
#' @param seed RNG seed
#' @param groups Optional group labels for the training rows; members of a
#'   group stay in the same CV fold and the same calibration half
#' @param levels Confidence levels for which the quantile is stored
#' @return List with \code{error_model}, \code{alpha_q} (at \code{confidence}),
#'   \code{confidence}, \code{levels} and \code{alpha_grid}
#' @export
train_normalized_conformal <- function(pca, response_model,
                                       X_train, y_train,
                                       confidence    = 0.68,
                                       committees    = NULL,
                                       k_cv          = 5L,
                                       conformal_subset_n = 10000L,
                                       seed          = 42L,
                                       groups        = NULL,
                                       levels        = .CONFORMAL_LEVELS) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required for normalized conformal calibration.")

  cm <- if (!is.null(committees)) committees
        else if (!is.null(response_model$committees)) response_model$committees
        else 5L

  set.seed(seed)
  T_full <- .project_pca(pca, X_train)
  n_full <- nrow(T_full)
  if (is.null(groups)) groups <- seq_len(n_full)

  # 1. Subsample: Safanelli reuses outer-loop CV predictions (free); a single
  #    split needs a fresh CV, so cap its cost.
  sub_idx <- if (n_full > conformal_subset_n) sort(sample.int(n_full, conformal_subset_n))
             else seq_len(n_full)
  T_in <- T_full[sub_idx, , drop = FALSE]
  y_in <- y_train[sub_idx]
  g_in <- groups[sub_idx]
  n_in <- nrow(T_in)

  # 2. Cross-validated residuals
  folds   <- random_kfold(n_in, k = k_cv, seed = seed, groups = g_in)
  yhat_cv <- numeric(n_in)
  for (te in folds) {
    tr <- setdiff(seq_len(n_in), te)
    mdl_f <- Cubist::cubist(x = T_in[tr, , drop = FALSE], y = y_in[tr],
                            committees = cm)
    yhat_cv[te] <- as.numeric(predict(mdl_f, T_in[te, , drop = FALSE]))
  }
  delta <- abs(y_in - yhat_cv)
  delta[!is.finite(delta)] <- NA_real_

  # 3. Error model on one half, scores on the other
  halves  <- .grouped_partition(n_in, c(0.5, 0.5), g_in, seed)
  fit_ix  <- halves[[1]][is.finite(delta[halves[[1]]])]
  cal_ix  <- halves[[2]]
  err_mdl <- Cubist::cubist(x = T_in[fit_ix, , drop = FALSE], y = delta[fit_ix],
                            committees = cm)

  delta_hat <- as.numeric(predict(err_mdl, T_in[cal_ix, , drop = FALSE]))
  delta_hat <- pmax(delta_hat, .Machine$double.eps)
  alpha     <- delta[cal_ix] / delta_hat
  alpha     <- alpha[is.finite(alpha)]

  alpha_grid <- .conformal_quantile(alpha, levels)

  list(
    error_model        = err_mdl,
    alpha_q            = .conformal_quantile(alpha, confidence),
    confidence         = confidence,
    levels             = levels,
    alpha_grid         = alpha_grid,
    n_calibration      = length(alpha),
    cv_folds           = k_cv,
    conformal_subset_n = n_in
  )
}

# Split-conformal quantile with the finite-sample correction: the
# ceiling((n + 1) q) / n empirical quantile guarantees coverage >= q.
.conformal_quantile <- function(alpha, q) {
  n   <- length(alpha)
  adj <- pmin(ceiling((n + 1) * q) / n, 1)
  as.numeric(stats::quantile(alpha, adj, type = 1, names = FALSE))
}

# The non-conformity quantile for a requested confidence, read from a model's
# stored grid. Older models stored a single level; they can only serve that one.
.conformal_alpha <- function(conf, confidence) {
  levels <- unlist(conf$levels)
  grid   <- unlist(conf$alpha_grid)
  if (length(levels) && length(grid) == length(levels)) {
    if (confidence < min(levels) - 1e-9 || confidence > max(levels) + 1e-9)
      stop(sprintf("confidence must be between %.2f and %.2f for this model.",
                   min(levels), max(levels)))
    return(stats::approx(levels, grid, xout = confidence, rule = 2)$y)
  }
  stored <- as.numeric(conf$confidence %||% 0.68)
  if (abs(confidence - stored) > 1e-9)
    warning(sprintf(paste0("This model was calibrated for %.0f%% intervals only; ",
                           "returning %.0f%% intervals, not %.0f%%. ",
                           "Retrain or download newer models for other levels."),
                    100 * stored, 100 * stored, 100 * confidence), call. = FALSE)
  as.numeric(conf$alpha_q)
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
