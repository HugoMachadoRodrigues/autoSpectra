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
  # Floor only real values. Uncovered bands stay NA: flooring them would turn a
  # missing band into absorbance -log(1e-6) = 13.8, a spike ten times the real
  # signal that then dominates SNV.
  eps <- 1e-6
  Mnum[which(is.finite(Mnum) & Mnum < eps)] <- eps
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
  keep <- mask == 1L & is.finite(M)
  nk   <- rowSums(keep)
  Mk   <- M; Mk[!keep] <- 0
  mu   <- rowSums(Mk) / pmax(nk, 1L)
  dev  <- (Mk - mu) * keep
  sg   <- sqrt(rowSums(dev^2) / pmax(nk, 1L))
  sg[!is.finite(sg) | sg == 0] <- 1
  out  <- dev / sg
  out[nk < 3L, ] <- 0          # too few covered bands to normalise
  dimnames(out) <- dimnames(M)
  out
}

# ---- Savitzky-Golay -------------------------------------------------------

#' Parse SG parameter string of the form "SG(w,p)", "SG(w,p,d)", etc.
#' @param st Character string like "SG_SMOOTH(11,2)": window length, polynomial
#'   order and (optionally) derivative order
#' @return Named list with m (window length), p, d
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
#' @param m Window length in bands (made odd if even)
#' @param p Polynomial order
#' @param ncols Number of spectral columns available
#' @return Adjusted list(m, p), or NULL if no valid window fits
sg_safe_params <- function(m, p, ncols) {
  if (!is.finite(p) || p < 0) p <- 2
  if (!is.finite(m) || m < 3) m <- 3
  m <- as.integer(m); p <- as.integer(p)
  if (m %% 2L == 0L) m <- m + 1L
  if (m <= p) m <- p + 1L + ((p + 1L) %% 2L == 0L)   # smallest odd window > p
  if (ncols < m) return(NULL)
  list(m = m, p = p)
}

# Least-squares Savitzky-Golay convolution weights for the centre point
.sg_weights <- function(w, p, d) {
  h <- (w - 1L) / 2L
  A <- outer(-h:h, 0:p, `^`)
  solve(crossprod(A), t(A))[d + 1L, ] * factorial(d)
}

# Smooth every row of a block that is fully finite. Interior points get the SG
# estimate; the (w-1)/2 points at each edge, where the window does not fit, keep
# their input value (smoothing) or the nearest interior value (derivative).
.sg_block <- function(B, wts, d) {
  B <- as.matrix(B)
  w <- length(wts); h <- (w - 1L) / 2L; n <- ncol(B)
  if (n < w) return(B)
  out <- B
  idx <- (h + 1L):(n - h)
  acc <- matrix(0, nrow(B), length(idx))
  for (j in seq_len(w)) acc <- acc + wts[j] * B[, idx - h - 1L + j, drop = FALSE]
  out[, idx] <- acc
  if (d > 0L) {
    out[, seq_len(h)]         <- acc[, 1L]
    out[, (n - h + 1L):n]     <- acc[, ncol(acc)]
  }
  out
}

#' Apply Savitzky-Golay filter row-wise to a spectral matrix
#'
#' Missing bands (\code{NA}) are respected: each contiguous run of finite bands
#' is filtered on its own and NAs stay NA, so an instrument that covers only
#' part of the range is smoothed exactly as it would be if measured alone.
#' Rows that share the same coverage are filtered together for speed.
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths)
#' @param m Window length in bands (odd; total window, not half-window)
#' @param p Polynomial order
#' @param d Derivative order (0 = smooth, 1 = first derivative)
#' @return Filtered matrix with the same dimensions and column names
#' @export
apply_sg_matrix <- function(M, m = 11, p = 2, d = 0) {
  M     <- as.matrix(M)
  pars  <- sg_safe_params(m, p, ncol(M))
  if (is.null(pars)) return(M)
  wts   <- .sg_weights(pars$m, pars$p, as.integer(d))
  fin   <- is.finite(M)
  out   <- M

  # Rows whose finite bands form one contiguous run [a, b] are the common case
  # (every OSSL instrument); filter them in blocks keyed by (a, b).
  first <- max.col(fin, ties.method = "first")
  last  <- ncol(M) + 1L - max.col(fin[, ncol(M):1, drop = FALSE], ties.method = "first")
  nfin  <- rowSums(fin)
  single_run <- nfin > 0L & nfin == (last - first + 1L)

  key <- ifelse(single_run, paste(first, last), NA_character_)
  for (k in unique(stats::na.omit(key))) {
    rows <- which(key == k)
    cols <- first[rows[1]]:last[rows[1]]
    out[rows, cols] <- .sg_block(M[rows, cols, drop = FALSE], wts, d)
  }

  # Rows with interior gaps: filter each finite run separately
  for (i in which(!single_run & nfin > 0L)) {
    r <- rle(fin[i, ])
    ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    for (j in which(r$values)) {
      cols <- starts[j]:ends[j]
      out[i, cols] <- .sg_block(M[i, cols, drop = FALSE], wts, d)
    }
  }
  out
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

# ---- Model input (shared by training and prediction) ---------------------

#' Turn resampled spectra into the matrix a model consumes
#'
#' The single place where training and prediction meet. Uncovered bands travel
#' through the whole pipeline as \code{NA} (absorbance, smoothing and SNV all
#' ignore them) and are set to 0 only at the very end. Filling them earlier is
#' what used to turn a missing band into a fake absorbance spike.
#'
#' @param X_res Matrix from \code{resample_to_grid()}; \code{NA} = not covered
#' @param steps Preprocessing steps (a family's \code{preprocess} field)
#' @param mask Optional coverage mask; derived from the \code{NA} pattern if
#'   omitted
#' @return List with \code{X} (model-ready, no \code{NA}) and \code{mask}
#' @export
prepare_model_input <- function(X_res, steps, mask = NULL) {
  X_res <- as.matrix(X_res)
  if (is.null(mask)) mask <- mask_from_coverage(X_res)
  X <- apply_pipeline(X_res, steps, mask = mask)
  X[!is.finite(X)] <- 0
  list(X = X, mask = mask)
}

# Low/high coverage limits of real instruments in OSSL v1.2, used both to
# simulate partial coverage in training and to score it afterwards.
.coverage_cuts <- function(sensor_type, wavegrid) {
  lo <- min(wavegrid); hi <- max(wavegrid)
  if (sensor_type == "visnir") {
    list(full_350_2500   = c(lo, hi),          # ASD
         from_400        = c(400, hi),         # Foss XDS, LUCAS
         from_800        = c(800, hi),         # Bruker MPA
         from_1000       = c(1000, hi),        # PerkinElmer
         from_1350       = c(1350, hi))        # NeoSpectra
  } else {
    list(full_600_4000   = c(lo, hi),
         from_650        = c(650, hi),
         from_700        = c(700, hi),
         to_3800         = c(max(lo, 600), 3800))
  }
}

#' Random spectral-dropout augmentation for masked-VAE training
#'
#' For every input row, generates copies that simulate partial-coverage
#' instruments by restricting the coverage mask to a randomly chosen instrument
#' range (VisNIR: from 350, 400, 800, 1000 or 1350 nm; MIR: 600, 650, 700 cm-1
#' or up to 3800 cm-1).
#'
#' Each copy has its uncovered bands removed and then goes through the full
#' preprocessing pipeline via \code{prepare_model_input()}, the same call
#' \code{predict_soil()} makes on a real spectrum. A simulated copy is
#' therefore exactly what the model would receive from an instrument with that
#' range. (Zeroing bands of an already-normalised spectrum, as earlier versions
#' did, gave inputs on a different scale from anything seen at prediction.)
#'
#' @param X Numeric matrix of resampled, NOT yet preprocessed spectra (rows =
#'   samples, cols = bands); uncovered bands \code{NA}
#' @param wavegrid Numeric vector of wavelength/wavenumber positions
#'   (length == ncol(X))
#' @param sensor_type "visnir" or "mir" (controls the instrument ranges)
#' @param n_copies How many augmented copies to generate per input row, in
#'   addition to the copy with the row's own coverage. Default 2.
#' @param seed Optional RNG seed for reproducibility
#' @param mask Optional true coverage mask of \code{X}; derived from its
#'   \code{NA} pattern if omitted. Augmented masks never exceed it.
#' @param steps Preprocessing pipeline (a family's \code{preprocess} field);
#'   default per-mask SNV only
#' @return List with \code{X} (model-ready copies, uncovered bands 0) and
#'   \code{mask} (same shape, 0/1)
#' @export
spectral_dropout <- function(X, wavegrid,
                             sensor_type = c("visnir", "mir"),
                             n_copies    = 2L,
                             seed        = NULL,
                             mask        = NULL,
                             steps       = "SNV_MASKED") {
  sensor_type <- match.arg(sensor_type)
  if (!is.null(seed)) set.seed(seed)

  X      <- as.matrix(X)
  if (is.null(mask)) mask <- mask_from_coverage(X)
  mask   <- as.matrix(mask)
  n_rows <- nrow(X)
  cuts   <- .coverage_cuts(sensor_type, wavegrid)
  in_cut <- lapply(cuts, function(cc) as.integer(wavegrid >= cc[1] & wavegrid <= cc[2]))

  X_out    <- vector("list", n_copies + 1L)
  mask_out <- vector("list", n_copies + 1L)
  X_out[[1]]    <- .model_input_under_mask(X, mask, steps)   # own coverage
  mask_out[[1]] <- mask

  for (k in seq_len(n_copies)) {
    pick <- sample.int(length(cuts), n_rows, replace = TRUE)
    mk   <- mask * do.call(rbind, in_cut[pick])
    X_out[[k + 1L]]    <- .model_input_under_mask(X, mk, steps)
    mask_out[[k + 1L]] <- mk
  }

  list(
    X    = do.call(rbind, X_out),
    mask = do.call(rbind, mask_out)
  )
}

# What predict_soil() would feed the model if the instrument had measured only
# the bands where mask == 1.
.model_input_under_mask <- function(X, mask, steps) {
  X[mask == 0L] <- NA_real_
  prepare_model_input(X, steps, mask = mask)$X
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
