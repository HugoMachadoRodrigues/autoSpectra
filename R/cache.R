# R/cache.R — In-memory model cache for autoSpectra
#
# Models are large keras objects. Loading them from disk on every prediction
# call is expensive. This module keeps loaded models in a package-level
# environment so they are fetched from disk exactly once per R session.
#
# Cache keys are "family_id/prop" (e.g. "OSSL_VisNIR/oc").

# Package-level environment — persists for the lifetime of the R session
.autoSpectra_cache <- new.env(parent = emptyenv())

#' Retrieve a soilVAE model from the in-memory cache
#'
#' On the first call for a given \code{family_id}/\code{prop} combination the
#' model is loaded from disk via \code{load_soilVAE()} and stored in the cache.
#' Subsequent calls return the cached object instantly.
#'
#' @param family_id Family identifier (\code{"OSSL_VisNIR"} or \code{"OSSL_MIR"})
#' @param prop Soil property name (e.g. \code{"oc"})
#' @param model_dir Root directory where trained models are stored
#' @return List with elements \code{model} (keras Model) and \code{scaler}
#'   (list with \code{mean} and \code{sd})
#' @export
get_cached_model <- function(family_id, prop,
                             model_dir = getOption("autoSpectra.model_dir",
                                                   "models")) {
  key <- paste0(family_id, "/", prop)
  if (!exists(key, envir = .autoSpectra_cache, inherits = FALSE)) {
    mdl_info <- load_soilVAE(family_id, prop, model_dir)
    assign(key, mdl_info, envir = .autoSpectra_cache)
  }
  get(key, envir = .autoSpectra_cache, inherits = FALSE)
}

#' Pre-load all soilVAE models for a family into memory
#'
#' Iterates over all (or selected) properties for \code{family_id} and loads
#' each model into the in-memory cache. Call this once at Shiny app startup to
#' eliminate per-prediction disk I/O.
#'
#' @param family_id Family identifier (\code{"OSSL_VisNIR"} or \code{"OSSL_MIR"})
#' @param model_dir Root directory where trained models are stored
#' @param properties Character vector of properties to pre-load.
#'   \code{NULL} (default) loads all properties defined for the family.
#' @param verbose Print progress messages (default \code{TRUE})
#' @return Invisible character vector of successfully loaded property names
#' @export
preload_ossl_models <- function(family_id,
                                model_dir  = getOption("autoSpectra.model_dir",
                                                       "models"),
                                properties = NULL,
                                verbose    = TRUE) {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("Package 'keras' is required for model preloading.")

  fam   <- get_family(family_id)
  props <- if (is.null(properties)) fam$properties else properties

  model_subdir <- file.path(model_dir, family_id, "models")
  loaded <- character(0)
  skipped <- character(0)

  if (verbose) message("autoSpectra: preloading [", fam$label, "] ...")

  for (prop in props) {
    h5_path <- file.path(model_subdir, paste0(prop, ".h5"))
    if (!file.exists(h5_path)) {
      skipped <- c(skipped, prop)
      next
    }
    key <- paste0(family_id, "/", prop)
    if (exists(key, envir = .autoSpectra_cache, inherits = FALSE)) {
      loaded <- c(loaded, prop)   # already cached
      next
    }
    ok <- tryCatch({
      mdl_info <- load_soilVAE(family_id, prop, model_dir)
      assign(key, mdl_info, envir = .autoSpectra_cache)
      TRUE
    }, error = function(e) {
      if (verbose) message("  [SKIP] ", prop, ": ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) loaded <- c(loaded, prop)
  }

  if (verbose) {
    message("  Loaded  : ", length(loaded), " models")
    if (length(skipped) > 0)
      message("  Missing : ", length(skipped), " (not yet trained)")
  }
  invisible(loaded)
}

#' List models currently held in the in-memory cache
#'
#' @return Character vector of cache keys (\code{"family_id/prop"})
#' @export
list_cached_models <- function() {
  ls(envir = .autoSpectra_cache)
}

#' Clear the in-memory model cache
#'
#' Releases all cached keras model objects. Models will be reloaded from disk
#' on the next prediction call.
#'
#' @return Invisible \code{NULL}
#' @export
clear_model_cache <- function() {
  rm(list = ls(envir = .autoSpectra_cache), envir = .autoSpectra_cache)
  invisible(NULL)
}
