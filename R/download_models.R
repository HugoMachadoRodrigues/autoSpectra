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

# Zenodo record that hosts the pre-trained model zips. Use the CONCEPT record
# ID: download_ossl_models() resolves it to the latest version through the
# Zenodo API, so publishing a new version needs no code change.
#
# The models need a deposit of their own (a dataset). The software record that
# Zenodo's GitHub integration creates for each release (concept 19004686)
# holds only the source archive and cannot carry extra files.
.MODELS_ZENODO_RECORD <- "19004686"

# ---- Publishing allowlist ------------------------------------------------

# Filenames that may be published to Zenodo. This is an ALLOWLIST on purpose:
# a denylist ships whatever it fails to recognise, which is exactly how the
# Cubist bundles came to carry the entire OSSL training set (see
# .strip_cubist_data / .strip_pca_scores). Every entry below is either model
# parameters or an aggregate statistic; nothing here is per-sample.
.PUBLISHABLE_PATTERNS <- c(
  "\\.keras$",                 # masked soilVAE weights
  "_cubist\\.rds$",            # Cubist response model + PCA basis (scores stripped)
  "_cubist_error\\.rds$",      # conformal error model (training blob stripped)
  "_scaler\\.rds$",            # target mean/sd
  "_scaler\\.json$",
  "_metrics\\.json$",          # per-property metrics, AD stats, conformal scalars
  "^metrics_summary\\.json$",  # per-family roll-up
  "^cv10_summary\\.json$",     # random 10-fold CV roll-up (per-fold metrics)
  "^lodo_summary\\.json$"      # leave-one-dataset-out roll-up (per-fold metrics)
)

#' Is this artefact cleared for publication?
#'
#' Applied by \code{upload_models.R} to every file before it is zipped for
#' Zenodo. Anything unmatched aborts the upload so a new artefact must be
#' reviewed deliberately rather than shipped by default.
#'
#' @param paths Character vector of file paths (only the basename is tested)
#' @return Logical vector, \code{TRUE} where the file may be published
#' @keywords internal
.is_publishable <- function(paths) {
  bn <- basename(paths)
  matched <- lapply(.PUBLISHABLE_PATTERNS, function(p) grepl(p, bn))
  Reduce(`|`, matched, init = logical(length(bn)))
}

# Files of the LATEST version of a Zenodo record (a concept ID resolves to
# its latest version), as a list of key / url pairs.
.zenodo_record_files <- function(record) {
  api  <- paste0("https://zenodo.org/api/records/", record)
  resp <- tryCatch(httr::GET(api, httr::timeout(60)),
                   error = function(e) stop("Could not reach Zenodo: ",
                                            conditionMessage(e), call. = FALSE))
  if (httr::http_error(resp))
    stop("Zenodo record ", record, " returned HTTP ", httr::status_code(resp),
         ".", call. = FALSE)
  rec   <- httr::content(resp, as = "parsed", simplifyVector = FALSE)
  files <- rec$files %||% list()
  list(
    id      = rec$id,
    version = rec$metadata$version %||% "?",
    keys    = vapply(files, function(f) f$key %||% "", character(1)),
    urls    = vapply(files, function(f) f$links$self %||% "", character(1))
  )
}

# Where to get <family>_models.zip from in a record. Either the zip is a file
# of the record itself (a manual deposit), or the record is the archive Zenodo's
# GitHub integration makes of the autoSpectra-models repository: a single zip
# whose top folder holds the family zips.
.locate_family_zip <- function(rec, zip_name) {
  hit <- match(zip_name, rec$keys)
  if (!is.na(hit)) return(list(url = rec$urls[hit], inner = NULL))
  zips <- which(grepl("\\.zip$", rec$keys))
  if (length(zips) == 1L) return(list(url = rec$urls[zips], inner = zip_name))
  stop("'", zip_name, "' is not in Zenodo record ", rec$id, " (version ",
       rec$version, "). It contains: ",
       if (length(rec$keys)) paste(rec$keys, collapse = ", ") else "no files",
       ".\nThe pre-trained models may not be published yet; to train them ",
       "yourself see train_ossl_models().", call. = FALSE)
}

# Pull <family>_models.zip out of a repository archive. Its path inside the
# archive starts with a folder named after the repository and commit.
.extract_inner_zip <- function(archive, zip_name, exdir) {
  members <- utils::unzip(archive, list = TRUE)$Name
  inner   <- members[basename(members) == zip_name]
  if (length(inner) != 1L)
    stop("'", zip_name, "' not found in the downloaded model archive.",
         call. = FALSE)
  utils::unzip(archive, files = inner, exdir = exdir, junkpaths = TRUE)
  file.path(exdir, zip_name)
}

#' Download the pre-trained OSSL models from Zenodo
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
#' @param family_id Character vector of family IDs to download: any of
#'   \code{"OSSL_VisNIR"}, \code{"OSSL_MIR"}, \code{"OSSL_VisNIR_Cubist"},
#'   \code{"OSSL_MIR_Cubist"}. Default: all four.
#' @param model_dir Local directory where models are stored.
#'   Default: \code{getOption("autoSpectra.model_dir", "models")}.
#' @param zenodo_record Zenodo record ID (string); a concept ID resolves to
#'   its latest version. Normally left at default.
#' @param overwrite Logical. If \code{TRUE}, re-download even if models are
#'   already present. Default \code{FALSE}.
#' @param timeout_sec Download timeout in seconds (default 3600 = 1 h).
#' @return Invisible character vector of directories written.
#' @export
download_ossl_models <- function(family_id   = c("OSSL_VisNIR", "OSSL_MIR",
                                                 "OSSL_VisNIR_Cubist",
                                                 "OSSL_MIR_Cubist"),
                                 model_dir   = getOption("autoSpectra.model_dir",
                                                         "models"),
                                 zenodo_record = .MODELS_ZENODO_RECORD,
                                 overwrite   = FALSE,
                                 timeout_sec = 3600L) {
  family_id <- match.arg(family_id,
                         choices  = c("OSSL_VisNIR", "OSSL_MIR",
                                      "OSSL_VisNIR_Cubist",
                                      "OSSL_MIR_Cubist"),
                         several.ok = TRUE)

  downloaded   <- character(0)
  rec          <- NULL      # Zenodo record metadata, fetched on first need
  archive_path <- NULL      # repository archive, downloaded at most once
  on.exit(if (!is.null(archive_path)) unlink(archive_path), add = TRUE)

  for (fid in family_id) {
    dest_dir <- file.path(model_dir, fid)

    # Check if models already present (count .keras, .h5, or _cubist.rds)
    model_subdir <- file.path(dest_dir, "models")
    artefacts <- if (dir.exists(model_subdir))
      c(list.files(model_subdir, pattern = "\\.(keras|h5)$"),
        list.files(model_subdir, pattern = "_cubist\\.rds$"))
    else character(0)

    if (!overwrite && length(artefacts) >= 1L) {
      message("autoSpectra: [", fid, "] already downloaded (",
              length(artefacts), " model artefacts in '", model_subdir, "').")
      message("  Use overwrite = TRUE to re-download.")
      downloaded <- c(downloaded, dest_dir)
      next
    }

    zip_name <- paste0(fid, "_models.zip")
    if (is.null(rec)) rec <- .zenodo_record_files(zenodo_record)
    src <- .locate_family_zip(rec, zip_name)

    message("autoSpectra: Downloading [", fid, "] from Zenodo record ",
            rec$id, " (version ", rec$version, ") ...")
    message("  Dest: ", dest_dir)

    # A repository archive holds every family: download it once per call
    if (!is.null(src$inner) && !is.null(archive_path)) {
      tmp_dl <- archive_path
    } else {
      tmp_dl <- tempfile(fileext = ".zip")
      resp <- tryCatch(
        httr::GET(src$url,
                  httr::write_disk(tmp_dl, overwrite = TRUE),
                  httr::progress(),
                  httr::timeout(timeout_sec)),
        error = function(e) {
          stop("Download failed: ", conditionMessage(e),
               "\nCheck your internet connection and that Zenodo record ",
               zenodo_record, " is accessible.", call. = FALSE)
        }
      )
      if (httr::http_error(resp))
        stop("Zenodo returned HTTP ", httr::status_code(resp), " for '",
             zip_name, "'.", call. = FALSE)
      if (!is.null(src$inner)) archive_path <- tmp_dl
    }
    tmp_zip <- if (is.null(src$inner)) tmp_dl
               else .extract_inner_zip(tmp_dl, zip_name, tempfile())

    # Extract zip into model_dir (preserves OSSL_VisNIR/ subfolder structure)
    message("  Extracting ...")
    dir_create(model_dir)
    utils::unzip(tmp_zip, exdir = model_dir)
    if (is.null(src$inner)) unlink(tmp_zip)

    n_extracted <- length(list.files(model_subdir,
                                     pattern = "\\.(keras|h5)$|_cubist\\.rds$",
                                     recursive = TRUE))
    message("  Done. ", n_extracted, " model files extracted to '",
            model_subdir, "'.")
    downloaded <- c(downloaded, dest_dir)
  }

  invisible(downloaded)
}

#' Check whether pre-trained models are available locally
#'
#' Counts both \code{.keras} (VAE), \code{.h5} (legacy VAE), and
#' \code{_cubist.rds} (Cubist baseline) artefacts.
#'
#' @param family_id Family ID
#' @param model_dir Root model directory
#' @return Named logical vector: \code{TRUE} if at least one trained model is found
#' @export
models_available <- function(family_id = c("OSSL_VisNIR", "OSSL_MIR",
                                           "OSSL_VisNIR_Cubist",
                                           "OSSL_MIR_Cubist"),
                             model_dir = getOption("autoSpectra.model_dir",
                                                   "models")) {
  family_id <- match.arg(family_id, several.ok = TRUE)
  vapply(family_id, function(fid) {
    d <- file.path(model_dir, fid, "models")
    if (!dir.exists(d)) return(FALSE)
    # At least one trained property. (A fixed minimum of 5 reported FALSE for
    # the published families, which ship 4 properties each.)
    n <- length(list.files(d, pattern = "\\.(keras|h5)$")) +
         length(list.files(d, pattern = "_cubist\\.rds$"))
    n >= 1L
  }, logical(1))
}
