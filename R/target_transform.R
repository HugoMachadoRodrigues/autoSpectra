# R/target_transform.R - log1p target handling with a bounded back-transform
#
# Skewed concentrations (oc, n.tot, ...) are modelled on a log1p scale. The
# naive inverse, expm1(), is exponential: a prediction that is wrong on the log
# scale becomes catastrophically wrong in original units.
#
# This is not hypothetical. Under leave-one-dataset-out validation, where the
# model is asked to predict a soil population it never saw, VisNIR Cubist scored
# (n-weighted R2 in original units):
#
#   n.tot    -14.83     log1p       one fold of 19,806 samples at R2 = -46.8
#   oc        -1.59     log1p       one fold of 19,747 samples at R2 =  -6.8
#   ph.h2o    +0.45     no log1p    holds up
#   clay.tot  +0.42     no log1p    holds up
#
# The split is exactly by transform, and the mechanism is expm1: on the same
# folds the model's ranking is still reasonable (median CCC 0.76-0.80), so the
# predictions are not random -- they are a few log-scale overshoots blown up
# into enormous absolute errors.
#
# The fix is to refuse to extrapolate. A model calibrated on organic carbon up
# to 78% has no basis for predicting 4,000%, and clamping to the observed
# training range on the LOG scale bounds the exponentiation at its source.
# Clamping after expm1() would be too late: the value has already exploded and
# the clamp would silently pile every runaway sample onto the ceiling.

#' Forward log1p transform for properties whose policy asks for it
#'
#' @param y Numeric vector of raw values
#' @param prop Property name; the transform is applied only when it appears in
#'   \code{ossl_l1_log_transform}
#' @return Numeric vector on the modelling scale, with values below -1 set NA
#' @export
apply_target_transform <- function(y, prop) {
  y <- suppressWarnings(as.numeric(y))
  if (!(prop %in% ossl_l1_log_transform)) return(y)
  y[!(is.finite(y) & y > -1)] <- NA_real_
  log1p(y)
}

#' Training-scale bounds used to keep a back-transform from extrapolating
#'
#' Recorded at fit time and stored beside the model, so prediction applies the
#' same bound the model was calibrated within.
#'
#' @param y_model Numeric vector of TRAINING targets on the modelling scale
#'   (already log1p-transformed where the policy applies)
#' @param pad Fractional slack on the observed range, so a sample marginally
#'   beyond the training extremes is not flattened (default 0.05 = 5% of range)
#' @return List with \code{lo} and \code{hi} on the modelling scale
#' @export
target_bounds <- function(y_model, pad = 0.05) {
  y <- y_model[is.finite(y_model)]
  if (!length(y)) return(list(lo = -Inf, hi = Inf))
  lo <- min(y); hi <- max(y)
  slack <- pad * (hi - lo)
  list(lo = lo - slack, hi = hi + slack)
}

#' Invert the target transform without letting it explode
#'
#' Clamps to the calibration range on the modelling scale BEFORE exponentiating,
#' then inverts. For untransformed properties the clamp still applies, which
#' costs nothing and keeps predictions physically plausible.
#'
#' @param y_hat Numeric vector of predictions on the modelling scale
#' @param prop Property name
#' @param bounds List from \code{target_bounds()}, or \code{NULL} to skip
#'   clamping (which restores the old, unbounded behaviour)
#' @param floor_zero Clamp the result at 0 for non-negative properties
#'   (default \code{TRUE}; concentrations cannot be negative)
#' @return Numeric vector in original units
#' @export
invert_target_transform <- function(y_hat, prop, bounds = NULL,
                                    floor_zero = TRUE) {
  y_hat <- as.numeric(y_hat)

  if (!is.null(bounds) && is.finite(bounds$lo) && is.finite(bounds$hi))
    y_hat <- pmin(pmax(y_hat, bounds$lo), bounds$hi)

  out <- if (prop %in% ossl_l1_log_transform) expm1(y_hat) else y_hat

  if (floor_zero && prop %in% ossl_l1_log_transform)
    out <- pmax(out, 0)

  out
}

#' Score predictions in original units, with the bounded back-transform
#'
#' The single place training code should compute metrics for a transformed
#' property, so every family scores the same way.
#'
#' @param y_true_model Observed values on the modelling scale
#' @param y_pred_model Predicted values on the modelling scale
#' @param prop Property name
#' @param bounds List from \code{target_bounds()} (typically built from the
#'   training targets), or \code{NULL}
#' @return List of metrics from \code{metrics_from_y()}, in original units
#' @export
metrics_original_scale <- function(y_true_model, y_pred_model, prop,
                                   bounds = NULL) {
  metrics_from_y(
    invert_target_transform(y_true_model, prop, bounds = NULL),
    invert_target_transform(y_pred_model, prop, bounds = bounds)
  )
}
