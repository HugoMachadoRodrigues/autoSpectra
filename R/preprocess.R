# R/preprocess.R - Spectral preprocessing pipeline for autoSpectra
#
# Pipeline steps:
#   ABSORBANCE      : convert reflectance to absorbance (-log(R))
#   SG_SMOOTH(m,p)  : Savitzky-Golay smoothing  (derivative order = 0)
#   SG_DERIV(m,p,d) : Savitzky-Golay derivative (derivative order = d, typ. 1)
#   SNV             : Standard Normal Variate (Safanelli 2025 method)
#   SNV_MASKED      : Per-mask SNV; mean/sd computed only on covered bands
#                     (used by masked-VAE families for partial-coverage instruments)
#
# Mask convention: a "mask" is a 0/1 matrix the same shape as the spectrum.
#   1 = band is covered by the source instrument (after resampling)
#   0 = band is missing / out of instrument coverage

# ---- Reflectance / Absorbance helpers ------------------------------------

#' Convert reflectance to absorbance
#'
#' Accepts reflectance in 0-1 range or as percent 0-100 (auto-rescaled).
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
#' Per Safanelli et al. 2025: SNV is the canonical OSSL preprocessing step,
#' selected through the SS4GG ring-trial experiment for its ability to remove
#' instrument-specific scaling and particle-size scattering effects.
#'
#' @param M Numeric matrix
#' @return SNV-transformed matrix
#' @export
apply_snv <- function(M) {
  M_c   <- sweep(M, 1, rowMeans(M, na.rm = TRUE), "-")
  denom <- sqrt(rowMeans(M_c^2, na.rm = TRUE))
  denom[denom == 0 | is.na(denom)] <- 1
  sweep(M_c, 1, denom, "/")
}

#' Apply SNV using only the masked (covered) bands of each row
#'
#' For each row, the mean and SD are computed using only the columns where
#' \code{mask == 1}. The same row is then z-normalised across all columns
#' (covered + uncovered), and uncovered positions are zero-filled so the
#' VAE encoder sees a numerically stable input.
#'
#' This is the SNV variant used by the masked-VAE families: it makes the
#' normalisation invariant to which instrument covered which bands.
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)
#' @param mask Binary matrix of the same shape (1 = covered, 0 = missing)
#' @return SNV matrix; uncovered positions set to 0
#' @export
apply_snv_masked <- function(M, mask) {
  if (!identical(dim(M), dim(mask)))
    stop("apply_snv_masked: M and mask must have identical dimensions.")
  M    <- as.matrix(M); mask <- as.matrix(mask)
  out  <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  for (i in seq_len(nrow(M))) {
    keep <- mask[i, ] == 1L & is.finite(M[i, ])
    if (sum(keep) < 3L) next
    vals     <- M[i, keep]
    mu       <- mean(vals)
    sg       <- sqrt(mean((vals - mu)^2))
    if (!is.finite(sg) || sg == 0) sg <- 1
    out[i, ] <- (M[i, ] - mu) / sg
    out[i, !keep] <- 0
  }
  colnames(out) <- colnames(M)
  out
}

# ---- Savitzky-Golay -------------------------------------------------------

#' Parse SG parameter string of the form "SG(m,p)", "SG(m,p,d)", etc.
#' @param st Character string like "SG(11,2,1)"
#' @return Named list with m, p, d
parse_sg <- function(st) {
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
  min_req <- p
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

# ---- Pipeline dispatcher --------------------------------------------------

#' Apply a spectral preprocessing pipeline
#'
#' Step strings:
#' \itemize{
#'   \item \code{"ABSORBANCE"}: reflectance -> absorbance (-ln R)
#'   \item \code{"SG_SMOOTH(m,p)"}: Savitzky-Golay smoothing (d=0)
#'   \item \code{"SG_DERIV(m,p,d)"}: Savitzky-Golay derivative (d>=1)
#'   \item \code{"SNV"}: Standard Normal Variate (Safanelli 2025)
#'   \item \code{"SNV_MASKED"}: per-mask SNV; requires \code{mask} argument
#' }
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)
#' @param steps Character vector of preprocessing step strings
#' @param mask Optional binary matrix for \code{SNV_MASKED} (same shape as M).
#'   If omitted while a \code{SNV_MASKED} step is present, a fully-covered
#'   mask is assumed.
#' @param absorbance_base10 Logical; use base-10 log for absorbance conversion
#' @return Preprocessed matrix
#' @export
apply_pipeline <- function(M, steps, mask = NULL, absorbance_base10 = FALSE) {
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
    } else if (identical(st, "SNV_MASKED")) {
      mk  <- if (is.null(mask)) matrix(1L, nrow(out), ncol(out)) else mask
      out <- apply_snv_masked(out, mk)
    }
  }
  out
}

# ---- Convenience wrappers (Safanelli pipeline) ---------------------------

#' Canonical Safanelli preprocessing for VisNIR (reflectance input)
#'
#' Per Safanelli 2025: SNV alone is sufficient and superior to derivatives.
#' The \code{ABSORBANCE} step is added for VisNIR because OSSL VisNIR data
#' is stored as reflectance.
#'
#' @param M Numeric matrix of reflectance values
#' @return SNV-preprocessed absorbance matrix
#' @export
preprocess_visnir <- function(M) {
  apply_pipeline(M, steps = c("ABSORBANCE", "SNV"))
}

#' Canonical Safanelli preprocessing for MIR (absorbance input)
#'
#' @param M Numeric matrix of absorbance values
#' @return SNV-preprocessed matrix
#' @export
preprocess_mir <- function(M) {
  apply_pipeline(M, steps = c("SNV"))
}

# ---- Mask helpers (for masked-VAE family) --------------------------------

#' Build a coverage mask from a resampled spectral matrix
#'
#' After \code{resample_to_grid()} samples that the source instrument did
#' not cover end up as \code{NA} (\code{rule = 1} in \code{stats::approx}).
#' This helper turns that \code{NA} pattern into a 0/1 mask.
#'
#' @param M_resampled Numeric matrix returned by \code{resample_to_grid()}
#' @return Integer matrix of the same shape (1 = covered, 0 = missing/NA)
#' @export
mask_from_coverage <- function(M_resampled) {
  mk <- matrix(1L, nrow = nrow(M_resampled), ncol = ncol(M_resampled))
  mk[!is.finite(M_resampled)] <- 0L
  mk
}

#' Random spectral-dropout augmentation for masked-VAE training
#'
#' For every input row, generates additional copies that simulate partial-
#' coverage instruments by zeroing out the mask outside randomised cut-off
#' wavelengths. Per OSSL v1.2, the realistic VisNIR cut-offs are 350, 400,
#' 800, 1000, 1350 nm (corresponding to ASD, Foss, Bruker MPA, PerkinElmer,
#' NeoSpectra); MIR cut-offs cluster around 600-700 cm-1.
#'
#' @param X Numeric matrix of preprocessed spectra (rows = samples, cols = bands)
#' @param wavegrid Numeric vector of wavelength/wavenumber positions
#'   (length == ncol(X))
#' @param sensor_type "visnir" or "mir" (controls realistic cut-off pool)
#' @param n_copies How many augmented copies to generate per input row
#'   (in addition to the original full-coverage copy). Default 2.
#' @param seed Optional RNG seed for reproducibility
#' @return List with \code{X} (rbind of all copies) and \code{mask}
#'   (same shape, 0/1)
#' @export
spectral_dropout <- function(X, wavegrid,
                             sensor_type = c("visnir", "mir"),
                             n_copies    = 2L,
                             seed        = NULL) {
  sensor_type <- match.arg(sensor_type)
  if (!is.null(seed)) set.seed(seed)

  X      <- as.matrix(X)
  n_rows <- nrow(X)
  n_col  <- ncol(X)

  # Realistic instrument low-end cut-offs found in OSSL v1.2
  cuts <- if (sensor_type == "visnir") {
    list(
      list(lo = min(wavegrid), hi = max(wavegrid)),         # full ASD
      list(lo = 400,           hi = max(wavegrid)),         # Foss XDS
      list(lo = 800,           hi = max(wavegrid)),         # Bruker MPA
      list(lo = 1000,          hi = max(wavegrid)),         # PerkinElmer
      list(lo = 1350,          hi = min(2550, max(wavegrid))) # NeoSpectra
    )
  } else {
    list(
      list(lo = min(wavegrid), hi = max(wavegrid)),
      list(lo = 650,           hi = max(wavegrid)),
      list(lo = 700,           hi = max(wavegrid)),
      list(lo = max(min(wavegrid), 600), hi = 3800)
    )
  }

  # Always include the full-coverage copy first
  X_out    <- vector("list", n_copies + 1L)
  mask_out <- vector("list", n_copies + 1L)
  X_out[[1]]    <- X
  mask_out[[1]] <- matrix(1L, n_rows, n_col)

  for (k in seq_len(n_copies)) {
    idx <- sample.int(length(cuts), n_rows, replace = TRUE)
    mk  <- matrix(0L, n_rows, n_col)
    for (i in seq_len(n_rows)) {
      cut <- cuts[[idx[i]]]
      mk[i, ] <- as.integer(wavegrid >= cut$lo & wavegrid <= cut$hi)
    }
    X_out[[k + 1L]]    <- X * mk      # zero out uncovered bands
    mask_out[[k + 1L]] <- mk
  }

  list(
    X    = do.call(rbind, X_out),
    mask = do.call(rbind, mask_out)
  )
}

#' Which rows carry any spectral signal at all
#'
#' The OSSL L0 tables have one row per soil LAYER, not per scan: a layer with no
#' scan on a given instrument still gets a row, with every spectral column
#' empty. That is 52.4% of the VisNIR table and 36.9% of MIR.
#'
#' Such rows MUST be dropped before training. They cannot be caught downstream
#' with \code{is.finite(rowSums(X))}, because the pipeline replaces non-finite
#' cells with 0 first: an all-NA row becomes an all-zero row, passes
#' \code{is.finite()}, and survives into training as a constant vector (constant
#' -1 after SNV). Every such row is then identical to every other, so identical
#' inputs force identical predictions and impose a hard, model-independent
#' ceiling on R2 -- which is why two unrelated algorithms previously landed
#' within 0.05 of each other on every property.
#'
#' Partial coverage is fine, and is precisely what the masked soilVAE exists to
#' handle; only rows with no coverage at all are unusable.
#'
#' A minimum-band threshold was considered and rejected on the data: coverage in
#' OSSL v1.2 is perfectly bimodal. Every non-empty row has at least 1,051 of
#' 1,076 VisNIR bands (the 1,051 floor is LUCAS/Metrohm starting at 400 nm
#' rather than 350) or all 1,701 MIR bands. Nothing sits in between, so
#' \code{> 0} and \code{>= 500} select exactly the same rows and a threshold
#' would be arbitrary complexity. Re-check this if a future OSSL release mixes
#' in instruments with genuinely sparse coverage.
#'
#' @param mask Integer/logical coverage matrix from \code{mask_from_coverage()}
#' @return Logical vector, \code{TRUE} for rows with at least one covered band
#' @export
has_spectral_coverage <- function(mask) {
  mask <- as.matrix(mask)
  rowSums(mask > 0L, na.rm = TRUE) > 0L
}
