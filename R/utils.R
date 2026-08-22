# R/utils.R — autoSpectra shared utility functions

# Internal helper: return y when x is NULL/empty
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

#' Extract wavelength/wavenumber positions from a data frame's column names
#'
#' Spectral columns are identified by having column names that are entirely
#' numeric (e.g. \code{"350"}, \code{"2500"}, \code{"600.5"}). Columns whose
#' names contain any non-numeric letter characters (e.g. \code{"oc"},
#' \code{"wr.33kPa"}, \code{"Soil_ID"}) are excluded.
#'
#' @param df A data frame with spectral columns named as numbers (nm or cm-1)
#' @param id_col Name of the sample identifier column to exclude
#' @return A list with \code{wl} (numeric positions) and \code{cols} (column names)
#' @export
get_wavelengths <- function(df, id_col = "Soil_ID") {
  cols <- setdiff(names(df), id_col)
  cols <- as.character(cols)
  # Parse the raw column name directly — only purely-numeric names (e.g. "350")
  # convert without NA.  Names like "wr.33kPa" or "clay.tot" stay as NA.
  wl  <- suppressWarnings(as.numeric(cols))
  idx <- !is.na(wl) & wl > 50   # spectral positions always > 50 nm / cm-1
  list(wl = wl[idx], cols = cols[idx])
}

#' Resample a spectral matrix to a target wavelength grid via linear interpolation
#'
#' Rows with fewer than 2 finite values (e.g., instruments that do not cover
#' the target range) are left as \code{NA} in the output and will be excluded
#' by the \code{is.finite(rowSums(X))} filter in the training pipeline.
#'
#' @param M Numeric matrix (rows = samples, cols = source wavelengths)
#' @param src_wl Numeric vector of source wavelength positions
#' @param target_wl Numeric vector of target wavelength positions
#' @return Numeric matrix resampled to target_wl
#' @export
resample_to_grid <- function(M, src_wl, target_wl) {
  M <- as.matrix(M)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = length(target_wl))
  for (i in seq_len(nrow(M))) {
    row_vals <- M[i, ]
    valid    <- is.finite(row_vals)
    if (sum(valid) < 2L) next    # leave row as NA; filter downstream
    out[i, ] <- stats::approx(
      x = src_wl[valid], y = row_vals[valid], xout = target_wl,
      rule = 1, ties = mean
    )$y
  }
  colnames(out) <- as.character(target_wl)
  out
}

#' Compute prediction metrics (RMSE, Bias, R², RPIQ, CCC)
#'
#' @param y_true Numeric vector of observed values
#' @param y_pred Numeric vector of predicted values
#' @return Named list: RMSE, Bias, R2, RPIQ, CCC (Lin's concordance)
#' @export
metrics_from_y <- function(y_true, y_pred) {
  ok     <- is.finite(y_true) & is.finite(y_pred)
  y_true <- y_true[ok]; y_pred <- y_pred[ok]
  rmse   <- sqrt(mean((y_true - y_pred)^2))
  bias   <- mean(y_pred - y_true)          # + = overestimate, − = underestimate
  sst    <- sum((y_true - mean(y_true))^2)
  sse    <- sum((y_true - y_pred)^2)
  r2     <- if (sst > 0) 1 - sse / sst else NA_real_
  rpiq   <- stats::IQR(y_true) / rmse
  # Lin's Concordance Correlation Coefficient (CCC)
  mu_t   <- mean(y_true);  mu_p   <- mean(y_pred)
  vt     <- mean((y_true - mu_t)^2); vp <- mean((y_pred - mu_p)^2)
  cov_tp <- mean((y_true - mu_t) * (y_pred - mu_p))
  denom  <- vt + vp + (mu_t - mu_p)^2
  ccc    <- if (denom > 0) 2 * cov_tp / denom else NA_real_
  list(RMSE = rmse, Bias = bias, R2 = r2, RPIQ = rpiq, CCC = ccc)
}

.normalize_target_transform <- function(method = c("none", "log1p", "sqrt", "auto")) {
  method <- if (is.null(method) || length(method) < 1L) "none" else method[1L]
  method <- tolower(as.character(method))
  match.arg(method, choices = c("none", "log1p", "sqrt", "auto"))
}

.ossl_log1p_target_properties <- c(
  "oc", "c.tot", "n.tot", "cec", "caco3",
  "al.ext", "al.ox", "al.dith",
  "ca.ext", "k.ext", "mg.ext", "na.ext", "p.ext",
  "fe.ext", "fe.dith", "fe.ox",
  "mn.ext", "zn.ext", "cu.ext", "b.ext",
  "s.tot", "s.ext", "acidity", "ec"
)

.ossl_none_target_properties <- c(
  "clay.tot", "silt.tot", "sand.tot",
  "ph.h2o", "ph.cacl2", "bd",
  "aggstb", "awc.33.1500kPa", "wr.33kPa", "wr.1500kPa"
)

.sample_skewness <- function(y) {
  y <- suppressWarnings(as.numeric(y))
  y <- y[is.finite(y)]
  if (length(y) < 3L)
    return(NA_real_)
  s <- stats::sd(y)
  if (!is.finite(s) || s <= 0)
    return(0)
  mean(((y - mean(y)) / s)^3)
}

resolve_target_transform <- function(prop = NULL, y = NULL,
                                     method = c("auto", "none", "log1p", "sqrt")) {
  method <- .normalize_target_transform(method)
  if (!identical(method, "auto"))
    return(method)

  prop <- if (is.null(prop) || length(prop) < 1L) NA_character_ else as.character(prop[1L])
  if (isTRUE(prop %in% .ossl_log1p_target_properties))
    return("log1p")
  if (isTRUE(prop %in% .ossl_none_target_properties))
    return("none")

  y <- suppressWarnings(as.numeric(y))
  y <- y[is.finite(y)]
  if (length(y) < 5L)
    return("none")
  if (any(y < 0))
    return("none")

  q <- suppressWarnings(stats::quantile(
    y,
    probs = c(0.25, 0.5, 0.75, 0.95, 0.99),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  q50 <- q[2L]
  q95 <- q[4L]
  q99 <- q[5L]
  skew <- .sample_skewness(y)
  zero_frac <- mean(y <= 0)
  ratio_99_50 <- if (is.finite(q50) && q50 > 0) q99 / q50 else Inf
  ratio_95_50 <- if (is.finite(q50) && q50 > 0) q95 / q50 else Inf

  if (zero_frac >= 0.05 && any(y > 0))
    return("log1p")
  if (is.finite(skew) && skew >= 2 && ratio_95_50 >= 4)
    return("log1p")
  if (is.finite(ratio_99_50) && ratio_99_50 >= 10)
    return("log1p")

  "none"
}

transform_target <- function(y, method = c("none", "log1p", "sqrt", "auto")) {
  method <- .normalize_target_transform(method)
  if (identical(method, "auto"))
    stop("Target transform 'auto' must be resolved before calling transform_target().")
  y <- suppressWarnings(as.numeric(y))
  if (identical(method, "none"))
    return(y)
  if (identical(method, "log1p")) {
    y[y < 0] <- NA_real_
    return(log1p(y))
  }
  y[y < 0] <- NA_real_
  sqrt(y)
}

inverse_target_transform <- function(y, method = c("none", "log1p", "sqrt", "auto")) {
  method <- .normalize_target_transform(method)
  if (identical(method, "auto"))
    stop("Target transform 'auto' must be resolved before calling inverse_target_transform().")
  y <- suppressWarnings(as.numeric(y))
  if (identical(method, "none"))
    return(y)
  if (identical(method, "log1p"))
    return(expm1(y))
  y ^ 2
}

make_target_bins <- function(y, n_bins = 4L) {
  y <- suppressWarnings(as.numeric(y))
  n_bins <- max(1L, suppressWarnings(as.integer(n_bins)[1L]))
  ok <- is.finite(y)
  out <- rep(NA_integer_, length(y))
  if (sum(ok) < 2L) {
    out[ok] <- 0L
    return(out)
  }
  probs <- seq(0, 1, length.out = n_bins + 1L)
  brks <- unique(as.numeric(stats::quantile(y[ok], probs = probs, na.rm = TRUE, type = 7)))
  if (length(brks) < 2L) {
    out[ok] <- 0L
    return(out)
  }
  bins <- cut(y[ok], breaks = brks, include.lowest = TRUE, labels = FALSE)
  bins[is.na(bins)] <- 1L
  out[ok] <- as.integer(bins - 1L)
  out
}

#' Split indices into train / calibration / test sets
#'
#' @param n Total number of samples
#' @param seed Random seed for reproducibility
#' @param p_train Fraction for training
#' @param p_cal Fraction for calibration (conformal prediction)
#' @param p_test Fraction for test
#' @return Named list with integer index vectors: train, calib, test
#' @export
split_idx <- function(n, seed = 42, p_train = 0.6, p_cal = 0.2, p_test = 0.2) {
  stopifnot(abs(p_train + p_cal + p_test - 1) < 1e-9)
  set.seed(seed)
  idx  <- sample.int(n)
  n_tr <- floor(p_train * n)
  n_ca <- floor(p_cal * n)
  tr   <- idx[seq_len(n_tr)]
  ca   <- idx[seq_len(n_ca) + n_tr]
  te   <- idx[-c(tr, ca)]
  list(train = tr, calib = ca, test = te)
}

#' Create a directory if it does not already exist
#' @param path Directory path to create
dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

#' Save a z-score scaler to RDS and JSON
#'
#' @param mean_y Mean of the target variable
#' @param sd_y Standard deviation of the target variable
#' @param path_base File path base (no extension)
save_scaler <- function(mean_y, sd_y, path_base, target_transform = "none") {
  scaler <- list(
    mean = as.numeric(mean_y),
    sd = as.numeric(sd_y),
    method = "z",
    target_transform = .normalize_target_transform(target_transform)
  )
  saveRDS(scaler, file = paste0(path_base, "_scaler.rds"))
  jsonlite::write_json(scaler, path = paste0(path_base, "_scaler.json"),
                       auto_unbox = TRUE, pretty = TRUE)
}

#' Compute per-band min/max feature range from a preprocessed matrix
#'
#' @param M Numeric matrix with named columns
#' @return Data frame with columns: band, min, max
feature_minmax <- function(M) {
  bands <- colnames(M)
  if (is.null(bands) || length(bands) != ncol(M)) {
    bands <- as.character(seq_len(ncol(M)))
  }
  data.frame(
    band = bands,
    min  = apply(M, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))),
    max  = apply(M, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))),
    check.names = FALSE
  )
}

#' Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles, latent stats)
#'
#' @param path_base File path base (no extension)
#' @param metrics List from metrics_from_y()
#' @param conf_q List with q90 and q95 conformal quantiles
#' @param feat_range Data frame from feature_minmax()
#' @param latent List with mu, Sigma, df, thr95
save_metrics_json <- function(path_base, metrics, conf_q, feat_range, latent) {
  obj <- list(
    metrics    = metrics,
    conformal  = list(q90 = conf_q[["q90"]], q95 = conf_q[["q95"]]),
    feat_range = feat_range,
    latent     = latent
  )
  jsonlite::write_json(obj, paste0(path_base, "_metrics.json"),
                       auto_unbox = TRUE, pretty = TRUE)
}

#' Safe shiny::validate wrapper
#' @param condition Logical condition to test
#' @param message Error message to show if condition is FALSE
vneed <- function(condition, message = "Validation failed.") {
  msg <- tryCatch({
    m <- paste(as.character(message), collapse = " ")
    if (!nzchar(m)) "Validation failed." else m
  }, error = function(e) "Validation failed.")
  shiny::validate(shiny::need(condition, msg))
}
