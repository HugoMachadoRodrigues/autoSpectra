# R/preprocess.R — Spectral preprocessing pipeline for autoSpectra
#
# Standard pipeline: ABSORBANCE -> SG_SMOOTH -> SG_DERIV
#   ABSORBANCE      : convert reflectance to absorbance (-log(R))
#   SG_SMOOTH(m,p)  : Savitzky-Golay smoothing   (derivative order = 0)
#   SG_DERIV(m,p,d) : Savitzky-Golay derivative  (derivative order = d, typically 1)
#   SNV             : Standard Normal Variate (kept for backward compat)
#   SG(m,p,d)       : Legacy single-step SG (backward compat)

#' Convert reflectance to absorbance
#'
#' Accepts reflectance in the range 0-1 or as percentage 0-100.
#' Values in the 0-100 range are automatically rescaled.
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths)
#' @param base10 Logical; if TRUE use -log10(R), otherwise -ln(R)
#' @return Absorbance matrix with same dimensions
#' @export
reflectance_to_absorbance <- function(M, base10 = FALSE) {
  Mnum <- as.matrix(M)
  if (is.finite(max(Mnum, na.rm = TRUE)) && max(Mnum, na.rm = TRUE) > 2)
    Mnum <- Mnum / 100
  eps <- 1e-6
  Mnum[Mnum < eps] <- eps
  if (base10) -log10(Mnum) else -log(Mnum)
}

#' Apply Standard Normal Variate (SNV) row-wise
#'
#' @param M Numeric matrix
#' @return SNV-transformed matrix
#' @export
apply_snv <- function(M) {
  M_c  <- sweep(M, 1, rowMeans(M, na.rm = TRUE), "-")
  denom <- sqrt(rowMeans(M_c^2, na.rm = TRUE))
  denom[denom == 0 | is.na(denom)] <- 1
  sweep(M_c, 1, denom, "/")
}

# ---- Internal SG helpers -----------------------------------------------

#' Parse SG parameter string of the form "SG(m,p)", "SG(m,p,d)", etc.
#' @param st Character string like "SG(11,2,1)"
#' @return Named list with m, p, d
parse_sg <- function(st) {
  # strip prefix up to and including first "("
  inner <- sub("^[^(]+\\(", "", st)
  inner <- sub("\\).*$", "", inner)
  parts <- trimws(strsplit(inner, ",", fixed = TRUE)[[1]])
  m <- suppressWarnings(as.integer(parts[1])); if (is.na(m)) m <- 11
  p <- suppressWarnings(as.integer(parts[2])); if (is.na(p)) p <- 2
  d <- if (length(parts) >= 3) suppressWarnings(as.integer(parts[3])) else 0
  if (is.na(d)) d <- 0
  list(m = m, p = p, d = d)
}

#' Compute safe SG parameters given available number of columns
#' @param m Half-window size (prospectr convention)
#' @param p Polynomial order
#' @param ncols Number of spectral columns
#' @return Adjusted list(m, p), or NULL if impossible
sg_safe_params <- function(m, p, ncols) {
  if (!is.finite(m) || m < 1) m <- 1
  if (!is.finite(p) || p < 0) p <- 2
  # window must be odd: 2m+1
  min_req <- p  # prospectr requires 2m+1 > p
  if (m < min_req) m <- min_req
  if (ncols < 3) return(NULL)
  max_m <- floor((ncols - 1) / 2)
  if (m > max_m) m <- max_m
  if (m < 1) return(NULL)
  list(m = m, p = p)
}

#' Apply Savitzky-Golay filter row-wise to a spectral matrix
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths)
#' @param m Half-window size (total window = 2m+1)
#' @param p Polynomial order
#' @param d Derivative order (0 = smooth, 1 = first derivative)
#' @return Filtered matrix; on error, returns input unchanged
#' @export
apply_sg_matrix <- function(M, m = 11, p = 2, d = 0) {
  ncols <- ncol(M)
  pars  <- sg_safe_params(m, p, ncols)
  if (is.null(pars)) return(M)
  m <- pars$m; p <- pars$p
  res <- t(apply(M, 1, function(r) {
    if (all(is.na(r))) return(r)
    tryCatch(
      prospectr::savitzkyGolay(r, m = m, p = p, w = d),
      error = function(e) r
    )
  }))
  colnames(res) <- colnames(M)
  res
}

# ---- Public pipeline function ------------------------------------------

#' Apply a spectral preprocessing pipeline
#'
#' The pipeline is defined as a character vector of step names:
#' \itemize{
#'   \item \code{"ABSORBANCE"}: convert reflectance → absorbance (−ln R)
#'   \item \code{"SG_SMOOTH(m,p)"}: Savitzky-Golay smooth (derivative = 0)
#'   \item \code{"SG_DERIV(m,p,d)"}: Savitzky-Golay derivative (d ≥ 1)
#'   \item \code{"SG(m,p,d)"}: legacy single-step SG (backward compatible)
#'   \item \code{"SNV"}: Standard Normal Variate
#' }
#'
#' The recommended two-step pipeline for OSSL VisNIR data is:
#' \code{c("ABSORBANCE", "SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)")}
#'
#' For MIR data (already in absorbance):
#' \code{c("SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)")}
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)
#' @param steps Character vector of preprocessing step strings
#' @param absorbance_base10 Logical; use base-10 log for absorbance conversion
#' @return Preprocessed matrix
#' @export
apply_pipeline <- function(M, steps, absorbance_base10 = FALSE) {
  out <- as.matrix(M)
  if (length(steps) == 0) return(out)
  for (st in steps) {
    if (identical(st, "ABSORBANCE")) {
      out <- reflectance_to_absorbance(out, base10 = absorbance_base10)
    } else if (startsWith(st, "SG_SMOOTH(")) {
      sg  <- parse_sg(st)
      out <- apply_sg_matrix(out, m = sg$m, p = sg$p, d = 0)
    } else if (startsWith(st, "SG_DERIV(")) {
      sg  <- parse_sg(st)
      out <- apply_sg_matrix(out, m = sg$m, p = sg$p, d = sg$d)
    } else if (startsWith(st, "SG(")) {
      sg  <- parse_sg(st)
      out <- apply_sg_matrix(out, m = sg$m, p = sg$p, d = sg$d)
    } else if (identical(st, "SNV")) {
      out <- apply_snv(out)
    }
  }
  out
}

#' Canonical two-step SG preprocessing for VisNIR (reflectance input)
#'
#' @param M Numeric matrix of reflectance values
#' @param m Half-window for SG (default 11)
#' @param p Polynomial order (default 2)
#' @return Preprocessed matrix: absorbance → smooth → 1st derivative
#' @export
preprocess_visnir <- function(M, m = 11, p = 2) {
  apply_pipeline(M, steps = c(
    "ABSORBANCE",
    sprintf("SG_SMOOTH(%d,%d)", m, p),
    sprintf("SG_DERIV(%d,%d,1)", m, p)
  ))
}

#' Canonical two-step SG preprocessing for MIR (absorbance input)
#'
#' @param M Numeric matrix of absorbance values (already in absorbance units)
#' @param m Half-window for SG (default 11)
#' @param p Polynomial order (default 2)
#' @return Preprocessed matrix: smooth → 1st derivative
#' @export
preprocess_mir <- function(M, m = 11, p = 2) {
  apply_pipeline(M, steps = c(
    sprintf("SG_SMOOTH(%d,%d)", m, p),
    sprintf("SG_DERIV(%d,%d,1)", m, p)
  ))
}
