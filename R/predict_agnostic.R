# R/predict_agnostic.R - sensor-agnostic prediction by instrument window
#
# A partial-coverage spectrum is predicted by the PCA-VAE family trained on
# the widest window it fully covers. Each window family is the OSSL PCA-VAE
# pipeline restricted to that window (R/registry.R); in the pilot they scored
# close to the full-coverage PCA-VAE on cut spectra, well above the masked
# soilVAE, which tried to cover every range with one network.

# Candidate families per sensor, widest window first.
.AGNOSTIC_WINDOWS <- list(
  visnir = c("OSSL_VisNIR_PCAVAE", "OSSL_VisNIR_PCAVAE_800",
             "OSSL_VisNIR_PCAVAE_1000", "OSSL_VisNIR_PCAVAE_1350"),
  mir    = c("OSSL_MIR_PCAVAE", "OSSL_MIR_PCAVAE_700_3800")
)

#' Choose the window family for each spectrum
#'
#' @param df Data frame as for \code{predict_soil()}
#' @param sensor \code{"visnir"} or \code{"mir"}
#' @return Character vector, one family id per row (\code{NA} when no window
#'   is fully covered)
#' @keywords internal
.agnostic_route <- function(df, sensor) {
  wl_info <- get_wavelengths(df, id_col = "Soil_ID")
  if (length(wl_info$wl) < 20)
    stop("No spectral columns detected (need at least 20 numeric columns).")
  X_src <- as.matrix(df[, wl_info$cols, drop = FALSE])
  route <- rep(NA_character_, nrow(df))
  for (fid in .AGNOSTIC_WINDOWS[[sensor]]) {
    todo <- is.na(route)
    if (!any(todo)) break
    X_res <- resample_to_grid(X_src[todo, , drop = FALSE], src_wl = wl_info$wl,
                              target_wl = model_registry[[fid]]$wavegrid)
    full <- rowSums(!is.finite(X_res)) == 0L
    route[which(todo)[full]] <- fid
  }
  route
}

#' Predict soil properties from spectra of any instrument range
#'
#' Sensor-agnostic prediction: each spectrum goes to the PCA-VAE model of the
#' widest instrument window it fully covers. VisNIR windows: 400-2500 nm (the
#' full-coverage PCA-VAE), 800-2500, 1000-2500 and 1350-2500 nm. MIR windows:
#' 600-4000 cm-1 (full coverage) and 700-3800 cm-1.
#'
#' @inheritParams predict_soil
#' @param sensor \code{"visnir"} or \code{"mir"}
#' @return Data frame as from \code{predict_soil()}, plus a
#'   \code{model_family} column naming the window model used for each row.
#'   Rows that cover no window get \code{NA} predictions and a warning.
#' @examples
#' \dontrun{
#' # Needs the window models: download_ossl_models(c("OSSL_VisNIR_PCAVAE",
#' #   "OSSL_VisNIR_PCAVAE_800", "OSSL_VisNIR_PCAVAE_1000", "OSSL_VisNIR_PCAVAE_1350"))
#' predict_soil_agnostic(spectra, sensor = "visnir")
#' }
#' @export
predict_soil_agnostic <- function(df, sensor = c("visnir", "mir"),
                                  properties = NULL,
                                  model_dir  = .default_model_dir(),
                                  with_pi    = FALSE,
                                  confidence = 0.68) {
  sensor <- match.arg(sensor)
  route  <- .agnostic_route(df, sensor)
  if (anyNA(route))
    warning(sprintf("%d spectra cover none of the %s windows (%s); left as NA.",
                    sum(is.na(route)), sensor,
                    paste(vapply(.AGNOSTIC_WINDOWS[[sensor]], function(f)
                      paste(range(model_registry[[f]]$wavegrid), collapse = "-"),
                      character(1)), collapse = ", ")),
            call. = FALSE)

  parts <- list()
  for (fid in unique(stats::na.omit(route))) {
    rows <- which(route == fid)
    p <- predict_soil(df[rows, , drop = FALSE], fid, properties = properties,
                      model_dir = model_dir, with_pi = with_pi,
                      confidence = confidence)
    p$.row <- rows
    parts[[fid]] <- p
  }
  out <- data.frame(Soil_ID = df[["Soil_ID"]], stringsAsFactors = FALSE)
  if (length(parts)) {
    all_p <- do.call(rbind, lapply(parts, function(p) p[, names(parts[[1]]), drop = FALSE]))
    for (nm in setdiff(names(all_p), c("Soil_ID", ".row"))) {
      out[[nm]] <- NA_real_
      out[[nm]][all_p$.row] <- all_p[[nm]]
    }
  }
  out$model_family <- route
  out
}
