# R/applicability.R - Applicability-domain methods for autoSpectra
#
# Two methods, one per algorithm family:
#
#   Cubist + PCA family : Q-statistic on PCA residuals
#                         (Jackson & Mudholkar 1979, used by Safanelli 2025)
#
#   VAE family          : Squared Mahalanobis distance in latent space
#                         (chi-squared threshold at df=latent_dim, alpha=0.05)
#
# Both produce a per-sample score and a TRUE/FALSE in-domain flag.

# ---- Q-statistic (PCA residuals) -----------------------------------------

#' Compute the Q-statistic for each row of a PCA residual matrix
#'
#' Q_i = e_i %*% e_i^T  (eq. 2 in Safanelli 2025).
#'
#' @param E Numeric matrix of PCA residuals (X - X_hat)
#' @return Numeric vector of Q values, length nrow(E)
#' @export
q_statistic <- function(E) {
  E <- as.matrix(E)
  rowSums(E * E)
}

#' Upper Confidence Limit for the Q-statistic (Jackson-Mudholkar 1979)
#'
#' Implements eq. 3 of Safanelli 2025:
#' \deqn{UCL_Q = \theta_1 \left[1 - \theta_2 h_0 \frac{1 - h_0}{\theta_1^2}
#'   + \frac{z_\alpha (2 \theta_2 h_0^2)^{1/2}}{\theta_1}\right]^{1/h_0}}
#' where \eqn{\theta_k = \mathrm{tr}(\mathrm{Cov}(E)^k)} and
#' \eqn{h_0 = 1 - \frac{2 \theta_1 \theta_3}{3 \theta_2^2}}.
#'
#' @param E Numeric matrix of training-set PCA residuals
#' @param alpha Significance level (default 0.01 = 99% confidence per Safanelli)
#' @return Numeric scalar UCL_Q
#' @export
q_statistic_ucl <- function(E, alpha = 0.01) {
  E    <- as.matrix(E)
  CovE <- stats::cov(E)
  ev   <- tryCatch(eigen(CovE, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) rep(0, ncol(CovE)))
  ev[ev < 0] <- 0
  th1 <- sum(ev)
  th2 <- sum(ev^2)
  th3 <- sum(ev^3)
  if (th1 <= 0 || th2 <= 0) return(NA_real_)
  h0  <- 1 - (2 * th1 * th3) / (3 * th2^2)
  if (!is.finite(h0) || h0 == 0) h0 <- 1
  z   <- stats::qnorm(1 - alpha)
  inner <- 1 - (th2 * h0 * (1 - h0)) / (th1^2) +
           (z * sqrt(2 * th2 * h0^2)) / th1
  if (!is.finite(inner) || inner <= 0) return(NA_real_)
  th1 * inner^(1 / h0)
}

#' Compute Q-statistic and in-domain flag for new spectra against a fitted PCA
#'
#' @param pca Fitted PCA list (from \code{.fit_pca()})
#' @param X_new Numeric matrix of new (preprocessed) spectra
#' @param ucl_q Upper Confidence Limit (from \code{q_statistic_ucl()})
#' @return Data frame with columns: q_stat, ucl_q, in_domain
#' @export
predict_q_statistic <- function(pca, X_new, ucl_q) {
  T_new   <- .project_pca(pca, X_new)
  X_hat   <- .reconstruct_pca(pca, T_new)
  E_new   <- as.matrix(X_new) - X_hat
  Q_new   <- q_statistic(E_new)
  data.frame(
    q_stat    = Q_new,
    ucl_q     = ucl_q,
    in_domain = Q_new <= ucl_q,
    stringsAsFactors = FALSE
  )
}

# ---- Mahalanobis distance in latent space (VAE family) -------------------

#' Squared Mahalanobis distance against a multivariate Gaussian (mu, Sigma)
#'
#' @param Z Numeric matrix of latent vectors (rows = samples, cols = latent dim)
#' @param mu Numeric vector (length = latent dim)
#' @param Sigma Numeric matrix (latent dim x latent dim)
#' @return Numeric vector of D^2 values, length nrow(Z)
#' @export
mahalanobis_d2 <- function(Z, mu, Sigma) {
  Z   <- as.matrix(Z)
  Sig_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  apply(Z, 1, function(z) {
    d <- z - mu
    as.numeric(t(d) %*% Sig_inv %*% d)
  })
}

#' Chi-squared threshold for the Mahalanobis test in latent space
#'
#' @param df Latent space dimension (degrees of freedom)
#' @param alpha Significance level (default 0.05)
#' @return Numeric scalar threshold
#' @export
mahalanobis_threshold <- function(df, alpha = 0.05) {
  stats::qchisq(1 - alpha, df = df)
}
