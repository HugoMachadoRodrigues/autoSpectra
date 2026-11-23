# R/download_models.R — Download pre-trained OSSL models from Zenodo
#
# Pre-trained models are hosted on Zenodo as zip archives, one per family
# (.PUBLISHED_FAMILIES): <family>_models.zip, or <family>_models.part<k>.zip
# when a family would pass GitHub's 100 MB file limit. Every zip holds
#   <family>/models/<prop>_...   (model parameters, scalers, metrics)
# and part 1 also holds the family-level summaries.
#
# Usage:
#   library(autoSpectra)
#   download_ossl_models()                         # default families
#   download_ossl_models("OSSL_MIR_PCAVAE")        # one family
#   download_ossl_models("all")                    # every published family
#   download_ossl_models(model_dir = "~/mymodels")

# Zenodo record that hosts the pre-trained model zips. Use the CONCEPT record
# ID: download_ossl_models() resolves it to the latest version through the
# Zenodo API, so publishing a new version needs no code change.
#
# The models live in github.com/HugoMachadoRodrigues/autoSpectra-models, which
# Zenodo archives as a dataset on every release (concept DOI
# 10.5281/zenodo.23950129). The package's own software record (concept
# 19004686) holds only the source code.
.MODELS_ZENODO_RECORD <- "23950129"

# GitHub repository whose releases carry each family zip as its own asset.
# The release tag equals the Zenodo version, so a family can be fetched alone
# (tens of MB) instead of through the full Zenodo archive (several GB).
.MODELS_GITHUB_REPO <- "HugoMachadoRodrigues/autoSpectra-models"

# ---- Publishing allowlist ------------------------------------------------

# Filenames that may be published to Zenodo. This is an ALLOWLIST on purpose:
# a denylist ships whatever it fails to recognise, which is exactly how the
# Cubist bundles came to carry the entire OSSL training set (see
# .strip_cubist_data / .strip_pca_scores). Every entry below is either model
# parameters or an aggregate statistic; nothing here is per-sample.
.PUBLISHABLE_PATTERNS <- c(
  "\\.keras$",                 # PCA-VAE ensemble weights
  "_pcavae\\.rds$",            # PCA-VAE: PCA basis (scores stripped) + member count
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

# Families published in the autoSpectra-models repository. The KSSL_*
# families are the OSSL Cubist / PCA-VAE pipelines trained on the KSSL rows
# only (Safanelli et al. 2025, "kssl" models). The *_PCAVAE_<window> families
# are the OSSL PCA-VAE restricted to one instrument window; predict_soil_agnostic()
# routes each spectrum to them. The masked soilVAE families (OSSL_VisNIR,
# OSSL_MIR) are no longer published: the window models are more accurate.
.PUBLISHED_FAMILIES <- c(
  "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist", "OSSL_NeoSpectra_Cubist",
  "OSSL_VisNIR_PCAVAE", "OSSL_MIR_PCAVAE", "OSSL_NeoSpectra_PCAVAE",
  "OSSL_VisNIR_PCAVAE_800", "OSSL_VisNIR_PCAVAE_1000",
  "OSSL_VisNIR_PCAVAE_1350", "OSSL_MIR_PCAVAE_700_3800",
  "KSSL_VisNIR_Cubist", "KSSL_MIR_Cubist",
  "KSSL_VisNIR_PCAVAE", "KSSL_MIR_PCAVAE"
)

# Downloaded when download_ossl_models() is called without family_id.
.DEFAULT_DOWNLOAD <- c("OSSL_VisNIR_PCAVAE", "OSSL_MIR_PCAVAE",
                       "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist")

# A family is one zip, <family>_models.zip, or, when it would pass GitHub's
# 100 MB file limit, several self-contained parts <family>_models.part<k>.zip.
# Returns the matching names in `keys` (basenames compared), parts in order.
.family_zip_names <- function(keys, family_id) {
  bn  <- basename(keys)
  pat <- paste0("^", gsub(".", "\\.", family_id, fixed = TRUE),
                "_models(\\.part([0-9]+))?\\.zip$")
  hit <- keys[grepl(pat, bn)]
  if (length(hit) < 2L) return(hit)
  k <- suppressWarnings(as.integer(sub(pat, "\\2", basename(hit))))
  hit[order(k, na.last = FALSE)]
}

# Split a family's files into zip parts of at most `max_bytes` (uncompressed),
# never separating the files of one property. `files` are paths under
# <family>/; family-level files (metrics_summary.json, ...) go in part 1.
# Returns a list of character vectors.
.partition_family_files <- function(files, max_bytes = 90e6) {
  bn    <- basename(files)
  props <- sub("_metrics\\.json$", "", bn[grepl("_metrics\\.json$", bn)])
  props <- props[order(nchar(props), decreasing = TRUE)]   # longest match first
  owner <- vapply(bn, function(b) {
    hit <- props[startsWith(b, paste0(props, "_")) | startsWith(b, paste0(props, "."))]
    if (length(hit)) hit[1] else ""
  }, character(1), USE.NAMES = FALSE)
  size  <- file.size(files)
  units <- split(seq_along(files), factor(owner, levels = unique(c("", sort(unique(owner))))))
  units <- units[lengths(units) > 0L]
  parts <- list(); cur <- integer(0); cur_size <- 0
  for (u in units) {
    us <- sum(size[u])
    if (length(cur) && cur_size + us > max_bytes) {
      parts[[length(parts) + 1L]] <- files[cur]; cur <- integer(0); cur_size <- 0
    }
    cur <- c(cur, u); cur_size <- cur_size + us
  }
  if (length(cur)) parts[[length(parts) + 1L]] <- files[cur]
  parts
}

# Where to get <family>_models.zip from in a record. Either the zip is a file
# of the record itself (a manual deposit), or the record is the archive Zenodo's
# GitHub integration makes of the autoSpectra-models repository: a single zip
# whose top folder holds the family zips.
.locate_family_zip <- function(rec, family_id) {
  direct <- .family_zip_names(rec$keys, family_id)
  if (length(direct))
    return(list(url = rec$urls[match(direct, rec$keys)], inner = NULL))
  zips <- which(grepl("\\.zip$", rec$keys))
  if (length(zips) == 1L) return(list(url = rec$urls[zips], inner = family_id))
  stop("'", family_id, "_models.zip' is not in Zenodo record ", rec$id, " (version ",
       rec$version, "). It contains: ",
       if (length(rec$keys)) paste(rec$keys, collapse = ", ") else "no files",
       ".\nThe pre-trained models may not be published yet; to train them ",
       "yourself see train_ossl_models().", call. = FALSE)
}

# Assets of the GitHub release `tag` of the models repository, with the
# SHA256 sums it publishes (SHA256SUMS asset). NULL when the release cannot
# be read (offline, rate limit, no such tag): the caller then falls back to
# the Zenodo archive.
.release_assets <- function(tag, repo = .MODELS_GITHUB_REPO) {
  tryCatch({
    api  <- sprintf("https://api.github.com/repos/%s/releases/tags/%s", repo, tag)
    resp <- httr::GET(api, httr::timeout(30),
                      httr::add_headers(Accept = "application/vnd.github+json"))
    if (httr::http_error(resp)) return(NULL)
    rel  <- httr::content(resp, as = "parsed", simplifyVector = FALSE)
    keys <- vapply(rel$assets, function(a) a$name %||% "", character(1))
    urls <- vapply(rel$assets, function(a) a$browser_download_url %||% "", character(1))
    if (!length(keys)) return(NULL)
    sums <- NULL
    i <- match("SHA256SUMS", keys)
    if (!is.na(i)) {
      r <- httr::GET(urls[i], httr::timeout(60))
      if (!httr::http_error(r))
        sums <- .parse_sha256sums(httr::content(r, as = "text", encoding = "UTF-8"))
    }
    list(keys = keys, urls = urls, sums = sums)
  }, error = function(e) NULL)
}

# "<hash>  <file>" lines -> named character vector (names = file basenames).
.parse_sha256sums <- function(txt) {
  ln <- strsplit(trimws(txt), "\r?\n")[[1]]
  ln <- ln[nzchar(ln)]
  m  <- regmatches(ln, regexec("^([0-9a-fA-F]{64})\\s+\\*?(.+)$", ln))
  m  <- m[lengths(m) == 3L]
  stats::setNames(tolower(vapply(m, `[`, "", 2L)),
                  basename(vapply(m, `[`, "", 3L)))
}

# SHA256 of a file as lowercase hex. Uses base R's tools::sha256sum()
# (R >= 4.5); openssl only on older R. (openssl's hashing segfaulted on Linux
# CI once Python's TensorFlow had loaded its own libssl in the same session.)
.file_sha256 <- function(path) {
  tools_ns <- asNamespace("tools")
  if (exists("sha256sum", envir = tools_ns, inherits = FALSE))
    return(tolower(unname(get("sha256sum", envir = tools_ns)(path))))
  con <- file(path, "rb"); on.exit(close(con))
  paste(as.character(unclass(openssl::sha256(con))), collapse = "")
}

# Stop unless `path` matches its published SHA256 (skipped when not listed).
.check_sha256 <- function(path, name, sums) {
  if (is.null(sums) || is.na(sums[name])) return(invisible(TRUE))
  got <- .file_sha256(path)
  if (!identical(tolower(got), unname(sums[name])))
    stop("Checksum mismatch for ", name, ": the download is corrupt. ",
         "Please retry.", call. = FALSE)
  invisible(TRUE)
}

# Download one URL to a local file, with progress. Kept separate so the tests
# can replace it and exercise download_ossl_models() without the network.
.download_to <- function(url, dest, timeout_sec, record) {
  resp <- tryCatch(
    httr::GET(url, httr::write_disk(dest, overwrite = TRUE),
              httr::progress(), httr::timeout(timeout_sec)),
    error = function(e) {
      stop("Download failed: ", conditionMessage(e),
           "\nCheck your internet connection and that Zenodo record ",
           record, " is accessible.", call. = FALSE)
    }
  )
  if (httr::http_error(resp))
    stop("Zenodo returned HTTP ", httr::status_code(resp), " for ", url, ".",
         call. = FALSE)
  invisible(dest)
}

# Pull a family's zip(s) out of a repository archive. Their path inside the
# archive starts with a folder named after the repository and commit.
.extract_inner_zip <- function(archive, family_id, exdir) {
  members <- utils::unzip(archive, list = TRUE)$Name
  inner   <- .family_zip_names(members, family_id)
  if (!length(inner))
    stop("'", family_id, "_models.zip' not found in the downloaded model archive.",
         call. = FALSE)
  utils::unzip(archive, files = inner, exdir = exdir, junkpaths = TRUE)
  file.path(exdir, basename(inner))
}

#' Download the pre-trained OSSL models from Zenodo
#'
#' Downloads the requested pre-trained model families and extracts them into
#' \code{model_dir}. The latest version is resolved through the Zenodo record;
#' each family is then fetched alone from the GitHub release of that version
#' (checked against its published SHA256 sums). If the release cannot be
#' reached, the full Zenodo archive (several GB) is downloaded instead. After the first
#' download, models are cached locally and loaded from disk (or from the
#' in-memory cache via \code{get_cached_model()}) on subsequent calls to
#' \code{predict_soil()}.
#'
#' @section Quick-start workflow:
#' \preformatted{
#' library(autoSpectra)
#' download_ossl_models()                  # one-time setup
#' predict_soil(df, "OSSL_VisNIR_PCAVAE")  # works immediately
#' }
#'
#' @param family_id Character vector of family IDs to download, any of:
#'   \code{"OSSL_VisNIR_Cubist"}, \code{"OSSL_MIR_Cubist"},
#'   \code{"OSSL_NeoSpectra_Cubist"}; \code{"OSSL_VisNIR_PCAVAE"},
#'   \code{"OSSL_MIR_PCAVAE"}, \code{"OSSL_NeoSpectra_PCAVAE"};
#'   \code{"OSSL_VisNIR_PCAVAE_800"}, \code{"OSSL_VisNIR_PCAVAE_1000"},
#'   \code{"OSSL_VisNIR_PCAVAE_1350"}, \code{"OSSL_MIR_PCAVAE_700_3800"}
#'   (instrument windows, used by \code{predict_soil_agnostic()});
#'   \code{"KSSL_VisNIR_Cubist"}, \code{"KSSL_MIR_Cubist"},
#'   \code{"KSSL_VisNIR_PCAVAE"}, \code{"KSSL_MIR_PCAVAE"} (trained on the KSSL
#'   rows only), or \code{"all"}. Default: the full-coverage VisNIR and MIR
#'   PCA-VAE and Cubist families. PCA-VAE families need the
#'   \pkg{soilVAE} package to predict.
#' @param model_dir Local directory where models are stored.
#'   Default: \code{getOption("autoSpectra.model_dir", "models")}.
#' @param zenodo_record Zenodo record ID (string); a concept ID resolves to
#'   its latest version. Normally left at default.
#' @param overwrite Logical. If \code{TRUE}, re-download even if models are
#'   already present. Default \code{FALSE}.
#' @param timeout_sec Download timeout in seconds (default 3600 = 1 h).
#' @return Invisible character vector of directories written.
#' @export
download_ossl_models <- function(family_id   = .DEFAULT_DOWNLOAD,
                                 model_dir   = getOption("autoSpectra.model_dir",
                                                         "models"),
                                 zenodo_record = .MODELS_ZENODO_RECORD,
                                 overwrite   = FALSE,
                                 timeout_sec = 3600L) {
  if (identical(family_id, "all")) family_id <- .PUBLISHED_FAMILIES
  family_id <- match.arg(family_id, choices = .PUBLISHED_FAMILIES, several.ok = TRUE)

  downloaded   <- character(0)
  rec          <- NULL      # Zenodo record metadata, fetched on first need
  rel          <- NULL      # GitHub release assets (FALSE when unavailable)
  archive_path <- NULL      # repository archive, downloaded at most once
  on.exit(if (!is.null(archive_path)) unlink(archive_path), add = TRUE)

  for (fid in family_id) {
    dest_dir <- file.path(model_dir, fid)

    # Check if models already present (count .keras, .h5, or _cubist.rds)
    model_subdir <- file.path(dest_dir, "models")
    artefacts <- if (dir.exists(model_subdir))
      c(list.files(model_subdir, pattern = "\\.(keras|h5)$"),
        list.files(model_subdir, pattern = "_(cubist|pcavae)\\.rds$"))
    else character(0)

    if (!overwrite && length(artefacts) >= 1L) {
      message("autoSpectra: [", fid, "] already downloaded (",
              length(artefacts), " model artefacts in '", model_subdir, "').")
      message("  Use overwrite = TRUE to re-download.")
      downloaded <- c(downloaded, dest_dir)
      next
    }

    if (is.null(rec)) rec <- .zenodo_record_files(zenodo_record)
    if (is.null(rel)) rel <- .release_assets(rec$version) %||% FALSE
    gh <- if (is.list(rel)) .family_zip_names(rel$keys, fid) else character(0)

    if (length(gh)) {
      # Only this family's zips, from the GitHub release of the same version
      message("autoSpectra: Downloading [", fid, "] (", length(gh),
              " file(s), models version ", rec$version, ") ...")
      message("  Dest: ", dest_dir)
      zips <- vapply(gh, function(k) {
        f <- tempfile(fileext = ".zip")
        .download_to(rel$urls[match(k, rel$keys)], f, timeout_sec, zenodo_record)
        .check_sha256(f, k, rel$sums)
        f
      }, character(1), USE.NAMES = FALSE)
      src <- NULL
    } else {
      src <- .locate_family_zip(rec, fid)
      message("autoSpectra: Downloading [", fid, "] from Zenodo record ",
              rec$id, " (version ", rec$version, ") ...")
      message("  Dest: ", dest_dir)
    }

    if (is.null(src)) {
      # already downloaded from the GitHub release
    } else if (is.null(src$inner)) {
      # Zips deposited as record files: one per family, or several parts
      zips <- vapply(src$url, function(u) {
        f <- tempfile(fileext = ".zip")
        .download_to(u, f, timeout_sec, zenodo_record)
        f
      }, character(1), USE.NAMES = FALSE)
    } else {
      # A repository archive holds every family: download it once per call
      if (is.null(archive_path)) {
        archive_path <- tempfile(fileext = ".zip")
        .download_to(src$url, archive_path, timeout_sec, zenodo_record)
      }
      zips <- .extract_inner_zip(archive_path, fid, tempfile())
    }

    # Extract into model_dir (each zip holds <family>/models/...)
    message("  Extracting ", length(zips), " archive(s) ...")
    dir_create(model_dir)
    for (z in zips) utils::unzip(z, exdir = model_dir)
    unlink(zips)

    n_extracted <- length(list.files(model_subdir,
                                     pattern = "\\.(keras|h5)$|_(cubist|pcavae)\\.rds$",
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
models_available <- function(family_id = .PUBLISHED_FAMILIES,
                             model_dir = getOption("autoSpectra.model_dir",
                                                   "models")) {
  family_id <- match.arg(family_id, several.ok = TRUE)
  vapply(family_id, function(fid) {
    d <- file.path(model_dir, fid, "models")
    if (!dir.exists(d)) return(FALSE)
    # At least one trained property. (A fixed minimum of 5 reported FALSE for
    # the published families, which ship 4 properties each.)
    n <- length(list.files(d, pattern = "\\.(keras|h5)$")) +
         length(list.files(d, pattern = "_(cubist|pcavae)\\.rds$"))
    n >= 1L
  }, logical(1))
}
