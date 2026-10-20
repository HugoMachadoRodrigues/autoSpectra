# R/predict.R - Prediction functions for autoSpectra
#
# Dispatches by family algorithm:
#   "vae_masked" -> masked-VAE (auto-builds mask from instrument coverage)
#   "vae"        -> standard VAE (full coverage required)
#   "cubist_pca" -> Cubist + PCA(120) (Safanelli replication)
#
# Predictions are automatically back-transformed for properties that were
# log1p-transformed at training time (per Safanelli 2025 Methods).

# ---- Top-level prediction ------------------------------------------------

#' Predict soil properties from a spectral data frame
#'
#' Resamples user spectra to the family's canonical grid, applies the family's
#' preprocessing pipeline (Safanelli's SNV for Cubist; SG-smooth + masked SNV
#' for masked VAE), generates a per-sample coverage mask from the resampled
#' \code{NA} pattern, and dispatches to the appropriate model loader.
#'
#' @param df Data frame with a \code{Soil_ID} column and numeric spectral
#'   columns named by wavelength (nm) or wavenumber (cm-1)
#' @param family_id Family identifier from \code{model_registry}
#' @param properties Properties to predict (default: all)
#' @param model_dir Root directory where models are stored
#' @param disable_pp Skip preprocessing (debug only)
#' @param with_pi Include prediction intervals when available (Cubist family)
#' @param confidence Confidence level for prediction intervals (default 0.68)
#' @return Data frame with \code{Soil_ID}, one column per predicted property,
#'   and (when \code{with_pi=TRUE}) \code{*_lo}/\code{*_hi} interval columns
#' @export
predict_soil <- function(df, family_id,
                         properties = NULL,
                         model_dir  = getOption("autoSpectra.model_dir",
                                                "models"),
                         disable_pp = FALSE,
                         with_pi    = FALSE,
                         confidence = 0.68) {

  fam   <- get_family(family_id)
  algo  <- fam$algorithm %||% "vae"
  props <- if (is.null(properties)) fam$properties else properties

  # Extract and resample spectra
  wl_info <- get_wavelengths(df, id_col = "Soil_ID")
  if (length(wl_info$wl) < 20)
    stop("No spectral columns detected (need at least 20 numeric columns).")

  X_src <- as.matrix(df[, wl_info$cols, drop = FALSE])
  X_res <- resample_to_grid(X_src, src_wl = wl_info$wl,
                            target_wl = fam$wavegrid)
  rownames(X_res) <- df[["Soil_ID"]]

  # Coverage mask from resampling NAs
  mask  <- mask_from_coverage(X_res)

  # Replace NA with 0 BEFORE preprocessing so SG/SNV doesn't propagate NAs
  X_res_nz <- X_res
  X_res_nz[!is.finite(X_res_nz)] <- 0

  X_proc <- if (disable_pp) X_res_nz
            else apply_pipeline(X_res_nz, fam$preprocess, mask = mask)

  out <- data.frame(Soil_ID = df[["Soil_ID"]], stringsAsFactors = FALSE)

  for (prop in props) {
    res <- tryCatch({
      if (identical(algo, "cubist_pca")) {
        .predict_one_cubist(family_id, prop, X_proc, model_dir,
                            with_pi = with_pi, confidence = confidence)
      } else if (identical(algo, "vae_masked")) {
        .predict_one_vae(family_id, prop, X_proc, model_dir, mask = mask)
      } else {
        .predict_one_vae(family_id, prop, X_proc, model_dir, mask = NULL)
      }
    }, error = function(e) list(y = rep(NA_real_, nrow(df))))

    out[[prop]] <- res$y
    if (with_pi && !is.null(res$pi_lo)) {
      out[[paste0(prop, "_lo")]] <- res$pi_lo
      out[[paste0(prop, "_hi")]] <- res$pi_hi
    }
  }

  # Back-transform to original units, bounded by the range each model was
  # calibrated on. An unbounded expm1() turns a log-scale overshoot on an
  # unfamiliar spectrum into an absurd concentration; see R/target_transform.R.
  for (prop in props) {
    bounds <- .load_target_bounds(family_id, prop, model_dir)
    out[[prop]] <- invert_target_transform(out[[prop]], prop, bounds = bounds)
    if (with_pi && !is.null(out[[paste0(prop, "_lo")]])) {
      out[[paste0(prop, "_lo")]] <-
        invert_target_transform(out[[paste0(prop, "_lo")]], prop, bounds = bounds)
      out[[paste0(prop, "_hi")]] <-
        invert_target_transform(out[[paste0(prop, "_hi")]], prop, bounds = bounds)
    }
  }
  out
}

#' Read a model's calibration range, if it was recorded at training time
#'
#' Models trained before the bounded back-transform have no stored range; they
#' return NULL and fall back to the old unbounded behaviour rather than failing.
#'
#' @keywords internal
.load_target_bounds <- function(family_id, prop, model_dir) {
  fp <- file.path(model_dir, family_id, "models", paste0(prop, "_scaler.rds"))
  if (!file.exists(fp)) return(NULL)
  sc <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(sc) || is.null(sc$bound_lo) || is.null(sc$bound_hi)) return(NULL)
  list(lo = as.numeric(sc$bound_lo), hi = as.numeric(sc$bound_hi))
}

# ---- Per-property dispatchers --------------------------------------------

.predict_one_vae <- function(family_id, prop, X_proc, model_dir, mask = NULL) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  base <- file.path(model_dir, family_id, "models", prop)
  has_keras <- file.exists(paste0(base, ".keras"))
  has_h5    <- file.exists(paste0(base, ".h5"))
  if (!has_keras && !has_h5)
    return(list(y = rep(NA_real_, nrow(X_proc))))

  mdl_info <- get_cached_model(family_id, prop, model_dir)
  mdl      <- mdl_info$model
  sc       <- mdl_info$scaler

  if (!is.null(mask)) {
    X_in <- .pack_masked_input(X_proc, mask)
  } else {
    X_in <- X_proc
  }

  expected <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]),
                       error = function(e) NA_integer_)
  if (!is.na(expected) && expected != ncol(X_in))
    stop("Input size mismatch for [", prop, "]: model expects ",
         expected, " columns, got ", ncol(X_in), ".")

  yhat_z <- .extract_prediction(mdl, X_in)
  mu  <- if (is.null(sc$mean)) 0 else sc$mean
  sg  <- if (is.null(sc$sd) || sc$sd == 0) 1 else sc$sd
  list(y = yhat_z * sg + mu)
}

.predict_one_cubist <- function(family_id, prop, X_proc, model_dir,
                                with_pi = FALSE, confidence = 0.68) {
  if (!requireNamespace("Cubist", quietly = TRUE))
    stop("Package 'Cubist' is required.")

  bundle_path <- file.path(model_dir, family_id, "models",
                           paste0(prop, "_cubist.rds"))
  if (!file.exists(bundle_path))
    return(list(y = rep(NA_real_, nrow(X_proc))))
  bundle <- readRDS(bundle_path)

  T_new  <- .project_pca(bundle$pca, X_proc)
  yhat   <- as.numeric(predict(bundle$model, T_new))

  res <- list(y = yhat)

  if (with_pi) {
    metrics_path <- file.path(model_dir, family_id, "models",
                              paste0(prop, "_metrics.json"))
    if (file.exists(metrics_path)) {
      info <- jsonlite::read_json(metrics_path)
      conf <- info$latent$conformal
      if (!is.null(conf) && !is.null(conf$alpha_q)) {
        # Reload error model from disk if present
        err_path <- file.path(model_dir, family_id, "models",
                              paste0(prop, "_cubist_error.rds"))
        err_mdl <- if (file.exists(err_path)) readRDS(err_path)
                   else conf$error_model
        if (!is.null(err_mdl)) {
          delta_hat <- as.numeric(predict(err_mdl, T_new))
          delta_hat <- pmax(delta_hat, 0)
          res$pi_lo <- yhat - delta_hat * as.numeric(conf$alpha_q)
          res$pi_hi <- yhat + delta_hat * as.numeric(conf$alpha_q)
        }
      }
    }
  }
  res
}

# ---- Applicability-domain dispatcher -------------------------------------

#' Compute applicability-domain scores for new spectra
#'
#' Dispatches to Q-statistic (Cubist family) or Mahalanobis on latent
#' space (VAE family) based on the family's algorithm.
#'
#' @param df Data frame with \code{Soil_ID} + spectral columns
#' @param family_id Family identifier
#' @param prop Soil property name (latent statistics are per-property for VAE)
#' @param model_dir Root model directory
#' @return Data frame with columns: \code{Soil_ID}, \code{score}, \code{threshold},
#'   \code{in_domain}, plus \code{method} attribute
#' @export
predict_applicability <- function(df, family_id, prop,
                                  model_dir = getOption("autoSpectra.model_dir",
                                                        "models")) {
  fam  <- get_family(family_id)
  algo <- fam$algorithm %||% "vae"

  wl_info <- get_wavelengths(df, id_col = "Soil_ID")
  X_src   <- as.matrix(df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, wl_info$wl, fam$wavegrid)
  mask    <- mask_from_coverage(X_res)
  X_res[!is.finite(X_res)] <- 0
  X_proc  <- apply_pipeline(X_res, fam$preprocess, mask = mask)

  if (identical(algo, "cubist_pca")) {
    bundle_path <- file.path(model_dir, family_id, "models",
                             paste0(prop, "_cubist.rds"))
    metrics_path <- file.path(model_dir, family_id, "models",
                              paste0(prop, "_metrics.json"))
    if (!file.exists(bundle_path) || !file.exists(metrics_path))
      stop("Missing Cubist bundle or metrics for ", family_id, "/", prop)
    bundle <- readRDS(bundle_path)
    info   <- jsonlite::read_json(metrics_path)
    ucl_q  <- as.numeric(info$latent$ucl_q)
    qd <- predict_q_statistic(bundle$pca, X_proc, ucl_q)
    out <- data.frame(
      Soil_ID    = df[["Soil_ID"]],
      score      = qd$q_stat,
      threshold  = qd$ucl_q,
      in_domain  = qd$in_domain,
      stringsAsFactors = FALSE
    )
    attr(out, "method") <- "q_statistic"
    return(out)
  }

  # VAE family - Mahalanobis on latent space
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  metrics_path <- file.path(model_dir, family_id, "models",
                            paste0(prop, "_metrics.json"))
  if (!file.exists(metrics_path))
    stop("Metrics file not found: ", metrics_path)

  info  <- jsonlite::read_json(metrics_path)
  mu_z  <- unlist(info$latent$mu)
  Sig_z <- do.call(rbind, lapply(info$latent$Sigma, unlist))
  thr95 <- info$latent$thr95 %||% mahalanobis_threshold(length(mu_z))

  mdl_info <- get_cached_model(family_id, prop, model_dir)
  mdl      <- mdl_info$model
  encoder  <- keras3::keras_model(inputs  = mdl$input,
                                  outputs = mdl$get_layer("z_mean")$output)

  X_in <- if (identical(algo, "vae_masked")) .pack_masked_input(X_proc, mask)
          else X_proc
  Z    <- predict(encoder, as.matrix(X_in), verbose = 0)
  d2   <- mahalanobis_d2(Z, mu_z, Sig_z)

  out <- data.frame(
    Soil_ID    = df[["Soil_ID"]],
    score      = d2,
    threshold  = thr95,
    in_domain  = d2 <= thr95,
    stringsAsFactors = FALSE
  )
  attr(out, "method") <- "mahalanobis"
  out
}

#' Round predictions and rename columns to display labels
#'
#' @param preds Data frame from \code{predict_soil()}
#' @param digits Number of decimal places (default 2)
#' @return Data frame with rounded values and human-readable column names
#' @export
format_predictions <- function(preds, digits = 2) {
  num_cols <- setdiff(names(preds), "Soil_ID")
  for (cl in num_cols)
    preds[[cl]] <- round(as.numeric(preds[[cl]]), digits)
  nice <- character(length(num_cols))
  for (i in seq_along(num_cols)) {
    nm <- num_cols[i]
    if (grepl("_lo$", nm)) {
      nice[i] <- paste0(property_label(sub("_lo$", "", nm)), " (lo)")
    } else if (grepl("_hi$", nm)) {
      nice[i] <- paste0(property_label(sub("_hi$", "", nm)), " (hi)")
    } else {
      nice[i] <- property_label(nm)
    }
  }
  names(preds) <- c("Soil_ID", nice)
  preds
}
