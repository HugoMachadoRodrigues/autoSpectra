# R/download_models.R — Download pre-trained OSSL soilVAE models from Zenodo
#
# Pre-trained models are hosted on Zenodo as zip archives, one per family:
#   OSSL_VisNIR_models.zip  (~300 MB)
#   OSSL_MIR_models.zip     (~500 MB)
#
# Each archive contains the full models/ sub-tree:
#   OSSL_VisNIR/
#     metrics_summary.json
#     models/
#       oc.h5  oc_scaler.rds  oc_metrics.json
#       clay.tot.h5  ...  (34 properties × 3 files)
#
# Usage:
#   library(autoSpectra)
#   download_ossl_models()                     # both families
#   download_ossl_models("OSSL_VisNIR")        # VisNIR only
#   download_ossl_models(model_dir = "~/mymodels")

# Zenodo record that hosts the pre-trained model zips.
# Updated each time new models are trained and uploaded.
.MODELS_ZENODO_RECORD <- "19004686"

#' Download pre-trained OSSL soilVAE models from Zenodo
#'
#' Downloads one or both pre-trained model archives from the autoSpectra
#' Zenodo deposit and extracts them into \code{model_dir}. After the first
#' download, models are cached locally and loaded from disk (or from the
#' in-memory cache via \code{get_cached_model()}) on subsequent calls to
#' \code{predict_soil()}.
#'
#' @section Quick-start workflow:
#' \preformatted{
#' library(autoSpectra)
#' download_ossl_models()                  # one-time setup
#' predict_soil(df, "OSSL_VisNIR")         # works immediately
#' }
#'
#' @param family_id Character vector of family IDs to download.
#'   One or both of \code{"OSSL_VisNIR"} and \code{"OSSL_MIR"}.
#'   Default: both families.
#' @param model_dir Local directory where models are stored.
#'   Default: \code{getOption("autoSpectra.model_dir", "models")}.
#' @param zenodo_record Zenodo record ID (string). Normally left at default.
#' @param overwrite Logical. If \code{TRUE}, re-download even if models are
#'   already present. Default \code{FALSE}.
#' @param timeout_sec Download timeout in seconds (default 3600 = 1 h).
#' @return Invisible character vector of directories written.
#' @export
download_ossl_models <- function(family_id   = c("OSSL_VisNIR", "OSSL_MIR"),
                                 model_dir   = getOption("autoSpectra.model_dir",
                                                         "models"),
                                 zenodo_record = .MODELS_ZENODO_RECORD,
                                 overwrite   = FALSE,
                                 timeout_sec = 3600L) {
  family_id <- match.arg(family_id,
                         choices  = c("OSSL_VisNIR", "OSSL_MIR"),
                         several.ok = TRUE)

  downloaded <- character(0)

  for (fid in family_id) {
    dest_dir <- file.path(model_dir, fid)

    # Check if models already present
    model_subdir <- file.path(dest_dir, "models")
    h5_files <- if (dir.exists(model_subdir))
      list.files(model_subdir, pattern = "\\.h5$")
    else character(0)

    if (!overwrite && length(h5_files) >= 5L) {
      message("autoSpectra: [", fid, "] already downloaded (",
              length(h5_files), " models in '", model_subdir, "').")
      message("  Use overwrite = TRUE to re-download.")
      downloaded <- c(downloaded, dest_dir)
      next
    }

    zip_name <- paste0(fid, "_models.zip")
    url      <- paste0("https://zenodo.org/record/", zenodo_record,
                       "/files/", zip_name, "?download=1")
    tmp_zip  <- tempfile(fileext = ".zip")

    message("autoSpectra: Downloading [", fid, "] from Zenodo record ",
            zenodo_record, " ...")
    message("  URL : ", url)
    message("  Dest: ", dest_dir)

    # Download with httr for progress reporting
    resp <- tryCatch(
      httr::GET(
        url,
        httr::write_disk(tmp_zip, overwrite = TRUE),
        httr::progress(),
        httr::timeout(timeout_sec)
      ),
      error = function(e) {
        stop("Download failed: ", conditionMessage(e),
             "\nCheck your internet connection and that Zenodo record ",
             zenodo_record, " is accessible.")
      }
    )

    if (httr::http_error(resp)) {
      status <- httr::status_code(resp)
      stop("Zenodo returned HTTP ", status, " for '", zip_name, "'.\n",
           "The models may not have been uploaded yet. ",
           "Train them locally with train_ossl.R first, then upload to Zenodo.")
    }

    # Extract zip into model_dir (preserves OSSL_VisNIR/ subfolder structure)
    message("  Extracting ...")
    dir_create(model_dir)
    utils::unzip(tmp_zip, exdir = model_dir)
    unlink(tmp_zip)

    n_extracted <- length(list.files(model_subdir, pattern = "\\.h5$",
                                     recursive = TRUE))
    message("  Done. ", n_extracted, " model files extracted to '",
            model_subdir, "'.")
    downloaded <- c(downloaded, dest_dir)
  }

  invisible(downloaded)
}

#' Check whether pre-trained models are available locally
#'
#' @param family_id Family ID (\code{"OSSL_VisNIR"} or \code{"OSSL_MIR"})
#' @param model_dir Root model directory
#' @return Named logical vector: \code{TRUE} if >= 5 .h5 files found
#' @export
models_available <- function(family_id = c("OSSL_VisNIR", "OSSL_MIR"),
                             model_dir = getOption("autoSpectra.model_dir",
                                                   "models")) {
  family_id <- match.arg(family_id, several.ok = TRUE)
  vapply(family_id, function(fid) {
    d <- file.path(model_dir, fid, "models")
    dir.exists(d) && length(list.files(d, pattern = "\\.h5$")) >= 5L
  }, logical(1))
}
