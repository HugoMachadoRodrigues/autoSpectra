# R/utils.R — autoSpectra shared utility functions

#' Null-coalescing operator: x %||% y returns y when x is NULL or empty.
#' @noRd
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

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
#' Target positions outside a row's measured range are \code{NA}, and so are
#' positions inside a gap in the measurements wider than \code{max_gap}: a
#' straight line drawn across a missing region is not a measurement, and the
#' coverage mask built from the result must not call it covered.
#'
#' Rows with fewer than 2 finite values are left entirely \code{NA}.
#'
#' @param M Numeric matrix (rows = samples, cols = source wavelengths)
#' @param src_wl Numeric vector of source wavelength positions
#' @param target_wl Numeric vector of target wavelength positions
#' @param max_gap Widest gap between two measured positions that is still
#'   interpolated across, in the units of \code{src_wl}. Default: 5 times the
#'   median source spacing.
#' @return Numeric matrix resampled to target_wl
#' @export
resample_to_grid <- function(M, src_wl, target_wl, max_gap = NULL) {
  M   <- as.matrix(M)
  ord <- order(src_wl)
  src_wl <- src_wl[ord]; M <- M[, ord, drop = FALSE]
  if (is.null(max_gap)) max_gap <- 5 * stats::median(diff(unique(src_wl)))
  out <- matrix(NA_real_, nrow = nrow(M), ncol = length(target_wl))
  for (i in seq_len(nrow(M))) {
    row_vals <- M[i, ]
    valid    <- is.finite(row_vals)
    if (sum(valid) < 2L) next    # leave row as NA; filter downstream
    xs <- src_wl[valid]
    y  <- stats::approx(x = xs, y = row_vals[valid], xout = target_wl,
                        rule = 1, ties = mean)$y
    gaps <- which(diff(xs) > max_gap)
    for (g in gaps) y[target_wl > xs[g] & target_wl < xs[g + 1L]] <- NA_real_
    out[i, ] <- y
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

#' Split indices into train / calibration / test sets
#'
#' With \code{groups}, every member of a group lands in the same set. Use it
#' when several rows are the same physical sample: in OSSL VisNIR, LUCAS
#' revisited the same points in 2009 and 2015, and a random split that puts one
#' visit in training and the other in test measures recall, not prediction.
#'
#' @param n Total number of samples
#' @param seed Random seed for reproducibility
#' @param p_train Fraction for training
#' @param p_cal Fraction for calibration (conformal prediction)
#' @param p_test Fraction for test
#' @param groups Optional vector of length \code{n}; rows sharing a value are
#'   kept together. Set sizes then match the fractions approximately.
#' @return Named list with integer index vectors: train, calib, test
#' @export
split_idx <- function(n, seed = 42, p_train = 0.6, p_cal = 0.2, p_test = 0.2,
                      groups = NULL) {
  stopifnot(abs(p_train + p_cal + p_test - 1) < 1e-9)
  parts <- .grouped_partition(n, c(p_train, p_cal, p_test), groups, seed)
  list(train = parts[[1]], calib = parts[[2]], test = parts[[3]])
}

# Partition 1..n into consecutive shares of a random order of GROUPS, so no
# group is split. Without groups every row is its own group.
.grouped_partition <- function(n, fractions, groups = NULL, seed = 42) {
  set.seed(seed)
  if (is.null(groups)) groups <- seq_len(n)
  stopifnot(length(groups) == n)
  g_order <- sample(unique(groups))
  g_rank  <- match(groups, g_order)
  cum     <- cumsum(tabulate(g_rank, nbins = length(g_order)))
  edges   <- c(0, cumsum(fractions)[-length(fractions)]) * n
  g_part  <- findInterval(cum - 0.5, edges)
  row_part <- g_part[g_rank]
  ord <- order(g_rank)
  lapply(seq_along(fractions), function(j) ord[row_part[ord] == j])
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
#' @param bounds Optional list from \code{target_bounds()}: the calibration range on
#'   the modelling scale, stored so prediction can bound the back-transform the
#'   same way training scored it. \code{NULL} stores no range.
save_scaler <- function(mean_y, sd_y, path_base, bounds = NULL) {
  scaler <- list(mean = as.numeric(mean_y), sd = as.numeric(sd_y), method = "z")
  # Calibration range on the modelling scale, so prediction can bound the
  # back-transform the same way training scored it. Without this a deployed
  # model extrapolates through expm1() and returns absurd concentrations for
  # spectra unlike anything it was trained on (see R/target_transform.R).
  if (!is.null(bounds)) {
    scaler$bound_lo <- as.numeric(bounds$lo)
    scaler$bound_hi <- as.numeric(bounds$hi)
  }
  saveRDS(scaler, file = paste0(path_base, "_scaler.rds"))
  jsonlite::write_json(scaler, path = paste0(path_base, "_scaler.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

#' Compute per-band min/max feature range from a preprocessed matrix
#'
#' Falls back to positional band ids when the matrix has no column names --
#' otherwise \code{colnames()} returns NULL and \code{data.frame()} fails with
#' "arguments imply differing number of rows: 0, 2".
#'
#' @param M Numeric matrix, ideally with named columns
#' @return Data frame with columns: band, min, max
feature_minmax <- function(M) {
  M     <- as.matrix(M)
  bands <- colnames(M)
  if (is.null(bands)) bands <- as.character(seq_len(ncol(M)))
  data.frame(
    band = bands,
    min  = apply(M, 2, function(x) suppressWarnings(min(x, na.rm = TRUE))),
    max  = apply(M, 2, function(x) suppressWarnings(max(x, na.rm = TRUE))),
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

#' Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles, latent stats)
#'
#' @param path_base File path base (no extension)
#' @param metrics List from metrics_from_y()
#' @param conf_q List with q90 and q95 conformal quantiles
#' @param latent List with mu, Sigma, df, thr95
#'
#' The file holds metrics and aggregate model statistics only. Per-band
#' training ranges (\code{feature_minmax()}) were stored here until 0.4.0;
#' nothing read them, and each value was one training sample's extreme, so they
#' are no longer published.
save_metrics_json <- function(path_base, metrics, conf_q, latent) {
  obj <- list(
    metrics    = metrics,
    conformal  = list(q90 = conf_q[["q90"]], q95 = conf_q[["q95"]]),
    latent     = latent
  )
  # digits = NA: full precision. jsonlite's default rounds to 4 decimals, which
  # stored the VisNIR Q-statistic limit (~1e-5) as 0 and flagged every
  # spectrum as out of domain.
  jsonlite::write_json(obj, paste0(path_base, "_metrics.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
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
