# R/utils.R — autoSpectra shared utility functions

#' Extract wavelength/wavenumber positions from a data frame's column names
#'
#' @param df A data frame with spectral columns named as numbers (nm or cm-1)
#' @param id_col Name of the sample identifier column to exclude
#' @return A list with `wl` (numeric positions) and `cols` (column names)
#' @export
get_wavelengths <- function(df, id_col = "Soil_ID") {
  cols <- setdiff(names(df), id_col)
  cols <- as.character(cols)
  wl   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", cols)))
  idx  <- !is.na(wl)
  list(wl = wl[idx], cols = cols[idx])
}

#' Resample a spectral matrix to a target wavelength grid via linear interpolation
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
    out[i, ] <- stats::approx(
      x = src_wl, y = M[i, ], xout = target_wl,
      rule = 1, ties = mean
    )$y
  }
  colnames(out) <- as.character(target_wl)
  out
}

#' Compute prediction metrics (RMSE, R², RPIQ)
#'
#' @param y_true Numeric vector of observed values
#' @param y_pred Numeric vector of predicted values
#' @return Named list with RMSE, R2, RPIQ
#' @export
metrics_from_y <- function(y_true, y_pred) {
  rmse <- sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
  sst  <- sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
  sse  <- sum((y_true - y_pred)^2, na.rm = TRUE)
  r2   <- if (sst > 0) 1 - sse / sst else NA_real_
  rpiq <- stats::IQR(y_true, na.rm = TRUE) / rmse
  list(RMSE = rmse, R2 = r2, RPIQ = rpiq)
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
save_scaler <- function(mean_y, sd_y, path_base) {
  scaler <- list(mean = as.numeric(mean_y), sd = as.numeric(sd_y), method = "z")
  saveRDS(scaler, file = paste0(path_base, "_scaler.rds"))
  jsonlite::write_json(scaler, path = paste0(path_base, "_scaler.json"),
                       auto_unbox = TRUE, pretty = TRUE)
}

#' Compute per-band min/max feature range from a preprocessed matrix
#'
#' @param M Numeric matrix with named columns
#' @return Data frame with columns: band, min, max
feature_minmax <- function(M) {
  data.frame(
    band = colnames(M),
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
