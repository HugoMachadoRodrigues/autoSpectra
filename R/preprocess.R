# R/preprocess.R — Spectral preprocessing pipeline for autoSpectra
#
# Standard pipeline: ABSORBANCE -> SG_SMOOTH -> SNV
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

#' Apply SNV using only valid (unmasked) spectral bands
#'
#' Computes mean and standard deviation exclusively over the valid bands
#' indicated by \code{mask}, then normalises those bands in place.
#' Missing bands (mask == 0) are set to zero in the output.
#' This is the correct SNV formulation for partial spectra (e.g. a NIR sensor
#' that only covers 1000–2500 nm on a 400–2500 nm grid).
#'
#' @param M    Numeric matrix (rows = samples, cols = wavelengths)
#' @param mask Integer or logical matrix, same shape as \code{M}.
#'   1 (or TRUE) = band is present; 0 (or FALSE) = band is absent.
#' @return SNV-normalised matrix; absent bands are set to 0.
#' @export
apply_snv_masked <- function(M, mask) {
  mask   <- as.matrix(mask + 0)            # coerce to numeric (0/1)
  M_v    <- ifelse(mask > 0, M, 0)         # zero out absent bands even when M has NA
  n_v    <- pmax(rowSums(mask), 1L)        # valid band count per row
  mu     <- rowSums(M_v) / n_v            # per-row mean (valid bands only)
  M_c    <- sweep(M_v, 1, mu, "-") * mask  # centred; absent bands stay 0
  sig    <- sqrt(rowSums(M_c^2) / pmax(n_v - 1L, 1L))  # unbiased SD
  sig[!is.finite(sig) | sig == 0] <- 1
  out    <- sweep(M_c, 1, sig, "/")
  if (!is.null(colnames(M))) colnames(out) <- colnames(M)
  out
}

.default_spectral_dropout_profiles <- function() {
  data.frame(
    profile_id = c("nir_1000_2500", "sw_nir_1350_2500"),
    label = c("NIR (1000-2500 nm)", "SW-NIR (1350-2500 nm)"),
    wl_lo = c(1000, 1350),
    wl_hi = c(2500, 2500),
    stringsAsFactors = FALSE
  )
}

.normalize_spectral_dropout_profiles <- function(profiles) {
  if (is.null(profiles)) {
    profiles <- .default_spectral_dropout_profiles()
  }
  profiles <- as.data.frame(profiles, stringsAsFactors = FALSE)
  if (!all(c("wl_lo", "wl_hi") %in% names(profiles))) {
    stop("profiles must contain columns 'wl_lo' and 'wl_hi'.")
  }
  if (!"profile_id" %in% names(profiles)) {
    profiles$profile_id <- paste0("profile_", seq_len(nrow(profiles)))
  }
  if (!"label" %in% names(profiles)) {
    profiles$label <- paste0(profiles$wl_lo, "-", profiles$wl_hi, " nm")
  }
  profiles$wl_lo <- as.numeric(profiles$wl_lo)
  profiles$wl_hi <- as.numeric(profiles$wl_hi)
  profiles
}

#' Sensor-agnostic spectral dropout via fixed profiles or random masking
#'
#' Creates augmented copies of the absorbance matrix using either a fixed bank
#' of physically plausible spectral profiles (`mode = "profiles"`, default) or
#' random lower-bound masks (`mode = "random"`). Bands outside the retained
#' range are zeroed out and SNV is applied only over the valid bands. One
#' additional full-spectrum copy is always included.
#'
#' @param X_abs     Absorbance or smoothed-absorbance matrix prior to SNV
#'   (rows = samples, cols = bands).
#' @param wavegrid  Numeric vector of wavelength positions matching \code{ncol(X_abs)}.
#' @param mode      Augmentation mode: \code{"profiles"} (default) or
#'   \code{"random"}.
#' @param profiles  Optional data frame of fixed spectral profiles with columns
#'   \code{wl_lo} and \code{wl_hi}. When omitted, uses the default masked-VAE
#'   bank: 1000-2500 nm and 1350-2500 nm.
#' @param n_aug     Number of random augmented copies to generate (default 9).
#'   Used only when \code{mode = "random"}.
#' @param wl_lo_min Minimum lower-bound cutoff to sample (default 400 nm).
#' @param wl_lo_max Maximum lower-bound cutoff to sample (default 1500 nm).
#'   Used only when \code{mode = "random"}.
#' @param seed      Optional RNG seed for reproducibility (default NULL).
#' @return Named list:
#'   \describe{
#'     \item{X}{Augmented SNV matrix}
#'     \item{mask}{Binary mask matrix, same shape as \code{X}}
#'     \item{cutoffs}{Numeric vector of lower-bound cutoffs for masked copies}
#'     \item{profiles}{Data frame describing the masked copies that were used}
#'     \item{mode}{Augmentation mode actually used}
#'     \item{n_copies}{Number of copies per sample, including full-spectrum}
#'   }
#' @export
spectral_dropout <- function(X_abs, wavegrid,
                             mode      = c("profiles", "random"),
                             profiles  = NULL,
                              n_aug     = 9L,
                              wl_lo_min = 400,
                              wl_lo_max = 1500,
                              seed      = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X_abs)
  d <- length(wavegrid)

  # Always include full-spectrum (no masking)
  mask_full <- matrix(1L, nrow = n, ncol = d)
  X_list    <- list(apply_snv_masked(X_abs, mask_full))
  mask_list <- list(mask_full)

  if (identical(mode, "profiles")) {
    profiles_used <- .normalize_spectral_dropout_profiles(profiles)
    for (i in seq_len(nrow(profiles_used))) {
      mk <- matrix(0L, nrow = n, ncol = d)
      mk[, wavegrid >= profiles_used$wl_lo[i] & wavegrid <= profiles_used$wl_hi[i]] <- 1L
      X_list    <- c(X_list,    list(apply_snv_masked(X_abs, mk)))
      mask_list <- c(mask_list, list(mk))
    }
    cutoffs <- profiles_used$wl_lo
  } else {
    # n_aug random lower-bound cutoffs, each applied to all samples
    cutoffs <- stats::runif(n_aug, min = wl_lo_min, max = wl_lo_max)
    for (co in cutoffs) {
      mk <- matrix(0L, nrow = n, ncol = d)
      mk[, wavegrid >= co] <- 1L
      X_list    <- c(X_list,    list(apply_snv_masked(X_abs, mk)))
      mask_list <- c(mask_list, list(mk))
    }
    profiles_used <- data.frame(
      profile_id = paste0("random_", seq_len(n_aug)),
      label = paste0(round(cutoffs), "-2500 nm"),
      wl_lo = cutoffs,
      wl_hi = max(wavegrid),
      stringsAsFactors = FALSE
    )
  }

  X_aug    <- do.call(rbind, X_list)
  mask_aug <- do.call(rbind, mask_list)
  rownames(X_aug) <- rownames(mask_aug) <- NULL
  if (!is.null(colnames(X_abs)))
    colnames(X_aug) <- colnames(mask_aug) <- colnames(X_abs)

  list(
    X = X_aug,
    mask = mask_aug,
    cutoffs = cutoffs,
    profiles = profiles_used,
    mode = mode,
    n_copies = length(X_list)
  )
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
#' The recommended VisNIR pipeline for multi-sensor OSSL data is:
#' \code{c("ABSORBANCE", "SG_SMOOTH(11,2)", "SNV")}
#'
#' For MIR data (already in absorbance):
#' \code{c("SG_SMOOTH(11,2)", "SNV")}
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

#' Canonical VisNIR preprocessing for multi-sensor reflectance input
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
    "SNV"
  ))
}

#' Canonical MIR preprocessing for multi-instrument absorbance input
#'
#' @param M Numeric matrix of absorbance values (already in absorbance units)
#' @param m Half-window for SG (default 11)
#' @param p Polynomial order (default 2)
#' @return Preprocessed matrix: smooth → 1st derivative
#' @export
preprocess_mir <- function(M, m = 11, p = 2) {
  apply_pipeline(M, steps = c(
    sprintf("SG_SMOOTH(%d,%d)", m, p),
    "SNV"
  ))
}
