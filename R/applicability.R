# R/applicability.R - Applicability-domain methods for autoSpectra
#
# Methods per algorithm family:
#
#   Cubist + PCA and    : Q-statistic on PCA residuals (Jackson & Mudholkar
#   PCA-VAE families      1979, used by Safanelli 2025) plus Hotelling T2 on
#                         the PC scores; a spectrum is in domain when it
#                         passes both
#
#   VAE family          : Squared Mahalanobis distance in latent space
#                         (chi-squared threshold at df=latent_dim, alpha=0.05)
#
# Each produces a per-sample score and a TRUE/FALSE in-domain flag.

# ---- Q-statistic (PCA residuals) -----------------------------------------

#' Compute the Q-statistic for each row of a PCA residual matrix
#'
#' Q_i = e_i %*% e_i^T  (eq. 2 in Safanelli 2025).
#'
#' @param E Numeric matrix of PCA residuals (X - X_hat)
#' @return Numeric vector of Q values, length nrow(E)
#' @examples
#' E <- matrix(rnorm(20), nrow = 5)
#' q_statistic(E)
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
#' @examples
#' E <- matrix(rnorm(200), nrow = 50)
#' q_statistic_ucl(E)
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
#' @examples
#' Z <- matrix(rnorm(20), nrow = 10)
#' mahalanobis_d2(Z, mu = colMeans(Z), Sigma = stats::cov(Z))
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
#' @examples
#' mahalanobis_threshold(df = 32)
#' @export
mahalanobis_threshold <- function(df, alpha = 0.05) {
  stats::qchisq(1 - alpha, df = df)
}

# ---- Hotelling T2 (PCA scores) -------------------------------------------
#
# The Q-statistic measures the distance of a spectrum from the PCA subspace
# (what the 120 components cannot reproduce). A spectrum can lie inside the
# subspace and still sit far from the training spectra along the components
# themselves: an unusual mixture of known spectral features. Hotelling T2 is
# that second distance, the Mahalanobis distance of the PC scores. The two are
# used together (Jackson 1991; Westerhuis et al. 2000).

#' Hotelling T2 of PCA scores
#'
#' \eqn{T^2_i = \sum_k t_{ik}^2 / \lambda_k}, where \eqn{\lambda_k} is the
#' variance of component \eqn{k} in the data the PCA was fitted on.
#'
#' @param scores Numeric matrix of PC scores (rows = samples)
#' @param sdev Standard deviation of each component (length \code{ncol(scores)})
#' @return Numeric vector of T2 values, length \code{nrow(scores)}
#' @examples
#' S <- matrix(rnorm(30), nrow = 10)
#' hotelling_t2(S, sdev = apply(S, 2, sd))
#' @export
hotelling_t2 <- function(scores, sdev) {
  scores <- as.matrix(scores)
  if (length(sdev) != ncol(scores))
    stop("sdev must have one value per column of scores")
  rowSums(sweep(scores^2, 2, sdev^2, "/"))
}

#' F-distribution limit for Hotelling T2 of a new observation
#'
#' \deqn{T^2_{lim} = \frac{k (n - 1)(n + 1)}{n (n - k)} F_{1-\alpha}(k, n - k)}
#' (Jackson 1991; Nomikos and MacGregor 1995). It assumes normally
#' distributed scores; the models record an empirical limit instead when they
#' are trained, and this one is used only for models trained without it.
#'
#' @param n Number of spectra the PCA was fitted on
#' @param k Number of components
#' @param alpha Significance level (default 0.01)
#' @return Numeric scalar limit
#' @examples
#' hotelling_t2_ucl(n = 50000, k = 120)
#' @export
hotelling_t2_ucl <- function(n, k, alpha = 0.01) {
  if (!is.finite(n) || n <= k) return(NA_real_)
  k * (n - 1) * (n + 1) / (n * (n - k)) * stats::qf(1 - alpha, k, n - k)
}

# Component standard deviations of a fitted PCA. Recorded as `sdev` by
# .fit_pca(). Older models keep only `var_explained`; for a scaled PCA (unit
# variance bands) the total variance is the number of bands, which recovers
# them. NULL when they cannot be recovered.
.pca_sdev <- function(pca) {
  if (!is.null(pca$sdev)) return(as.numeric(pca$sdev))
  scaled <- !is.null(pca$scale) && !all(pca$scale == 1)
  if (scaled && !is.null(pca$var_explained))
    return(sqrt(as.numeric(pca$var_explained) * nrow(pca$rotation)))
  NULL
}

#' Compute Hotelling T2 and in-domain flag for new spectra against a fitted PCA
#'
#' @param pca Fitted PCA list (from \code{.fit_pca()})
#' @param X_new Numeric matrix of new (preprocessed) spectra
#' @param ucl_t2 Limit (empirical, from the model's metrics, or from
#'   \code{hotelling_t2_ucl()})
#' @return Data frame with columns \code{t2}, \code{ucl_t2}, \code{in_domain};
#'   \code{t2} is \code{NA} when the PCA has no component variances
#' @export
predict_hotelling_t2 <- function(pca, X_new, ucl_t2) {
  sdev <- .pca_sdev(pca)
  t2 <- if (is.null(sdev)) rep(NA_real_, nrow(as.matrix(X_new)))
        else hotelling_t2(.project_pca(pca, X_new), sdev)
  data.frame(t2 = t2, ucl_t2 = ucl_t2, in_domain = t2 <= ucl_t2,
             stringsAsFactors = FALSE)
}
