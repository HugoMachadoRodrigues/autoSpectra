# R/predict.R — Prediction functions for autoSpectra

#' Predict soil properties from a spectral data frame using saved soilVAE models
#'
#' @param df Data frame with a Soil_ID column and spectral columns (numeric names
#'   in nm or cm-1)
#' @param family_id Family identifier (e.g., "OSSL_VisNIR", "ASD_DRY")
#' @param properties Character vector of properties to predict. Default: all
#'   available in the family.
#' @param model_dir Root directory where trained models are stored (default "models")
#' @param disable_pp Logical; skip preprocessing (for debug only)
#' @return Data frame with Soil_ID and one column per predicted property
#' @export
predict_soil <- function(df, family_id,
                         properties = NULL,
                         model_dir  = "models",
                         disable_pp = FALSE) {

  if (!requireNamespace("keras", quietly = TRUE))
    stop("Package 'keras' is required for prediction.")

  fam   <- get_family(family_id)
  props <- if (is.null(properties)) fam$properties else properties

  # Extract and resample spectra
  wl_info <- get_wavelengths(df, id_col = "Soil_ID")
  if (length(wl_info$wl) < 20)
    stop("No spectral columns detected (need at least 20 numeric columns).")

  X_src <- as.matrix(df[, wl_info$cols, drop = FALSE])
  X_res <- resample_to_grid(X_src, src_wl = wl_info$wl,
                            target_wl = fam$wavegrid)
  rownames(X_res) <- df[["Soil_ID"]]

  X_proc <- if (disable_pp) X_res else apply_pipeline(X_res, fam$preprocess)

  out <- data.frame(Soil_ID = df[["Soil_ID"]], stringsAsFactors = FALSE)

  for (prop in props) {
    fp_h5 <- file.path(model_dir, family_id, "models", paste0(prop, ".h5"))
    if (!file.exists(fp_h5)) {
      out[[prop]] <- NA_real_
      next
    }

    mdl_info <- load_soilVAE(family_id, prop, model_dir)
    mdl      <- mdl_info$model
    sc       <- mdl_info$scaler

    # Input shape check
    expected <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]),
                         error = function(e) NA_integer_)
    if (!is.na(expected) && expected != ncol(X_proc))
      stop("Input size mismatch for [", prop, "]: model expects ",
           expected, " bands, got ", ncol(X_proc),
           ". Check that the selected family matches the trained model.")

    yhat_z    <- .extract_prediction(mdl, X_proc)
    mu        <- if (is.null(sc$mean)) 0 else sc$mean
    sg        <- if (is.null(sc$sd) || sc$sd == 0) 1 else sc$sd
    out[[prop]] <- yhat_z * sg + mu
  }
  out
}

#' Compute applicability domain score for new samples
#'
#' Uses the Mahalanobis distance in latent space, calibrated against the
#' chi-squared threshold at 95% (saved during training).
#'
#' @param df Data frame with Soil_ID and spectral columns
#' @param family_id Family identifier
#' @param prop Soil property name (model must have been trained with latent stats)
#' @param model_dir Root model directory
#' @return Data frame with Soil_ID, mahal_dist, thr95, in_domain (logical)
#' @export
predict_applicability <- function(df, family_id, prop, model_dir = "models") {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("Package 'keras' is required.")

  fam <- get_family(family_id)

  # Load latent stats from metrics JSON
  base        <- file.path(model_dir, family_id, "models", prop)
  metrics_path <- paste0(base, "_metrics.json")
  if (!file.exists(metrics_path))
    stop("Metrics file not found: ", metrics_path)

  info   <- jsonlite::read_json(metrics_path)
  mu_z   <- unlist(info$latent$mu)
  Sig_z  <- do.call(rbind, lapply(info$latent$Sigma, unlist))
  thr95  <- info$latent$thr95

  # Get latent embeddings for new samples
  mdl_info <- load_soilVAE(family_id, prop, model_dir)
  mdl      <- mdl_info$model
  encoder  <- keras::keras_model(inputs = mdl$input,
                                 outputs = mdl$get_layer("latent")$output)

  wl_info <- get_wavelengths(df, id_col = "Soil_ID")
  X_src   <- as.matrix(df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, wl_info$wl, fam$wavegrid)
  X_proc  <- apply_pipeline(X_res, fam$preprocess)

  Z <- predict(encoder, as.matrix(X_proc), verbose = 0)

  # Mahalanobis distance
  Sig_inv <- tryCatch(solve(Sig_z), error = function(e) MASS::ginv(Sig_z))
  dists   <- apply(Z, 1, function(z) {
    d <- z - mu_z
    as.numeric(t(d) %*% Sig_inv %*% d)
  })

  data.frame(
    Soil_ID   = df[["Soil_ID"]],
    mahal_dist = dists,
    thr95      = thr95,
    in_domain  = dists <= thr95,
    stringsAsFactors = FALSE
  )
}

#' Round predictions and rename columns to fancy display labels
#'
#' @param preds Data frame from predict_soil()
#' @param digits Number of decimal places (default 2)
#' @return Data frame with rounded values and pretty column names
#' @export
format_predictions <- function(preds, digits = 2) {
  num_cols <- setdiff(names(preds), "Soil_ID")
  for (cl in num_cols)
    preds[[cl]] <- round(as.numeric(preds[[cl]]), digits)
  nice <- c("Soil_ID", vapply(num_cols, property_label, character(1)))
  names(preds) <- nice
  preds
}
