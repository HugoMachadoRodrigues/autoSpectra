# R/training_inputs.R - From joined OSSL data to what each trainer consumes
#
# One function builds the training inputs for every family, using the same
# preprocessing call as predict_soil(). Keeping this in the package (rather
# than in train_ossl.R) is what lets the tests check train/predict parity.

#' Group label per OSSL row: rows that are the same physical sample
#'
#' LUCAS revisited the same points in 2009 and 2015; OSSL keeps both visits as
#' separate rows whose \code{id.scan_local_c} is \code{"<year>.<point id>"}.
#' Those visits share a group. Every other row is its own group.
#'
#' @param joined_df Data frame from \code{ossl_prepare()}
#' @return Character vector, one label per row
#' @keywords internal
.ossl_sample_groups <- function(joined_df) {
  ids <- as.character(joined_df[["Soil_ID"]])
  ds  <- joined_df[["dataset.code_ascii_txt"]]
  loc <- joined_df[["id.scan_local_c"]]
  if (is.null(ds) || is.null(loc)) return(ids)
  lucas <- !is.na(ds) & ds == "LUCAS.SSL" & !is.na(loc) &
           grepl("^[0-9]{4}\\.", loc)
  ids[lucas] <- paste0("LUCAS:", sub("^[0-9]{4}\\.", "", loc[lucas]))
  ids
}

#' Build the training inputs for one model family
#'
#' @param fam Family entry from \code{model_registry}
#' @param joined_df Data frame from \code{ossl_prepare()}
#' @return List with
#'   \item{X}{model-ready matrix for Cubist families (\code{prepare_model_input()})}
#'   \item{X_raw}{for masked-VAE families: the resampled, unpreprocessed
#'     spectra with uncovered bands \code{NA} (what \code{train_maskedVAE()}
#'     expects); \code{NULL} otherwise}
#'   \item{mask}{coverage mask}
#'   \item{has_spectra}{rows with at least one covered band}
#'   \item{groups}{sample group per row (\code{.ossl_sample_groups()})}
#' @keywords internal
.family_training_inputs <- function(fam, joined_df) {
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)
  rm(X_src)
  mask <- mask_from_coverage(X_res)

  # OSSL L0 has one row per soil LAYER, not per scan, so ~half the rows carry
  # no spectrum at all; see has_spectral_coverage().
  has_spectra <- has_spectral_coverage(mask)
  groups      <- .ossl_sample_groups(joined_df)

  algo <- fam$algorithm %||% "vae"
  if (identical(algo, "vae_masked")) {
    # The masked VAE preprocesses inside the trainer, once per simulated
    # coverage, so it takes the resampled spectra as they are.
    return(list(X = NULL, X_raw = X_res, mask = mask,
                has_spectra = has_spectra, groups = groups))
  }
  mi <- prepare_model_input(X_res, fam$preprocess, mask = mask)
  list(X = mi$X, X_raw = NULL, mask = mask,
       has_spectra = has_spectra, groups = groups)
}
