# upload_models.R — Package trained models for Zenodo upload
#
# Run AFTER train_ossl.R completes. Creates per-family zip archives for
# the autoSpectra-models repository, which Zenodo archives on every release
# (concept DOI 10.5281/zenodo.23950129).
#
# Output, one per family in .PUBLISHED_FAMILIES (R/download_models.R):
#   dist/<family>_models.zip                 (family below 90 MB), or
#   dist/<family>_models.part<k>.zip         (larger families: GitHub refuses
#                                             files over 100 MB; every part is
#                                             a self-contained zip, and the
#                                             files of one property never span
#                                             two parts)
#
# Each zip contains:
#   <family_id>/models/<prop>.keras / <prop>_vae<k>.keras / <prop>_pcavae.rds
#   <family_id>/models/<prop>_cubist.rds / <prop>_cubist_error.rds
#   <family_id>/models/<prop>_scaler.rds / <prop>_metrics.json
#   <family_id>/metrics_summary.json (part 1)
#
# Every file passes the publishing allowlist, and every zip is audited for
# per-sample OSSL data; an archive that fails is deleted.
#
# Publishing procedure:
#   1. Copy the zip archives into a clone of
#      https://github.com/HugoMachadoRodrigues/autoSpectra-models and refresh
#      SHA256SUMS (shasum -a 256 *.zip > SHA256SUMS)
#   2. Commit, push, and create a GitHub release (e.g. v0.5.0)
#   3. Zenodo archives the release as a new version of the same concept record;
#      download_ossl_models() follows the latest version automatically
#
# Usage:
#   Rscript upload_models.R                              # every published family
#   Rscript upload_models.R OSSL_VisNIR                  # one family
#   Rscript upload_models.R --model-dir models_release   # models trained elsewhere

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION")) {
    if (requireNamespace("pkgload", quietly = TRUE))
      pkgload::load_all(".", quiet = TRUE)
    else
      for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
  } else {
    library(autoSpectra)
  }
})

if (!requireNamespace("zip", quietly = TRUE))
  stop("Package 'zip' is required: install.packages('zip')")

# The publishing allowlist lives in the package (R/download_models.R) so it is
# unit-tested. It is internal, so reach it through ::: when autoSpectra was
# attached rather than loaded from source.
is_publishable <- if (exists(".is_publishable", inherits = TRUE))
  get(".is_publishable", inherits = TRUE) else
  utils::getFromNamespace(".is_publishable", "autoSpectra")

pkg_internal <- function(name) if (exists(name, inherits = TRUE))
  get(name, inherits = TRUE) else utils::getFromNamespace(name, "autoSpectra")
partition_family_files <- pkg_internal(".partition_family_files")

args      <- commandArgs(trailingOnly = TRUE)
MODEL_DIR <- "models"
if (length(i <- which(args == "--model-dir"))) {
  MODEL_DIR <- args[i[1] + 1L]; args <- args[-c(i[1], i[1] + 1L)]
}
FAMILIES  <- if (length(args) > 0) args else pkg_internal(".PUBLISHED_FAMILIES")
DIST_DIR  <- "dist"
MAX_PART  <- 90e6    # bytes before compression; GitHub refuses files > 100 MB

dir_create(DIST_DIR)

audit_failures <- character(0)

for (fid in FAMILIES) {
  src_dir <- file.path(MODEL_DIR, fid)

  if (!dir.exists(src_dir)) {
    message("SKIP [", fid, "]: '", src_dir, "' not found (run train_ossl.R first).")
    next
  }

  # Enumerate model artefacts
  model_subdir <- file.path(src_dir, "models")
  keras_files  <- list.files(model_subdir, pattern = "\\.keras$",  full.names = FALSE)
  rds_files    <- list.files(model_subdir, pattern = "\\.rds$",    full.names = FALSE)
  json_files   <- list.files(model_subdir, pattern = "\\.json$",   full.names = FALSE)
  summary_json <- file.exists(file.path(src_dir, "metrics_summary.json"))

  n_models <- max(length(keras_files), length(rds_files))
  message("[", fid, "]: ", n_models, " models  |  ",
          length(keras_files), " .keras  |  ", length(rds_files), " .rds  |  ",
          length(json_files), " .json  |  metrics_summary: ", summary_json)

  if (n_models == 0) {
    message("  No model files found — skipping.")
    next
  }

  # Collect the files to publish. `.is_publishable()` is an ALLOWLIST: a
  # denylist ships whatever it does not recognise, which is precisely how the
  # Cubist bundles came to carry the whole OSSL training set — nothing in the
  # publisher objected. Unmatched files abort the upload instead.
  all_files <- list.files(src_dir, recursive = TRUE, full.names = TRUE)
  keep      <- is_publishable(all_files)

  if (any(!keep)) {
    message("  REFUSING to publish [", fid, "]: unrecognised file(s) under ",
            src_dir, ":")
    for (f in all_files[!keep])
      message("    - ", sub(paste0("^", src_dir, "/"), "", f),
              "  (", round(file.size(f) / 1e6, 2), " MB)")
    message("  Every published file must be model parameters or an aggregate ",
            "statistic. If this artefact is safe to distribute, add a pattern ",
            "to .PUBLISHABLE_PATTERNS in R/download_models.R; if it holds ",
            "per-sample OSSL data, do not ship it.")
    next
  }
  all_files <- all_files[keep]

  # A directory copied from another family (e.g. models_kssl/OSSL_MIR_Cubist
  # staged as KSSL_MIR_Cubist) still names the old family in its summary.
  summ <- file.path(src_dir, "metrics_summary.json")
  if (file.exists(summ)) {
    sid <- jsonlite::read_json(summ)$family_id
    if (!is.null(sid) && !identical(sid, fid)) {
      message("  REFUSING to publish [", fid, "]: metrics_summary.json says '",
              sid, "'. Set its family_id to '", fid, "' when staging the copy.")
      next
    }
  }

  # Remove this family's archives from an earlier run first, so a family that
  # now needs fewer parts does not leave stale ones behind.
  old_zips <- list.files(DIST_DIR, full.names = TRUE,
                         pattern = paste0("^", gsub(".", "\\.", fid, fixed = TRUE),
                                          "_models(\\.part[0-9]+)?\\.zip$"))
  unlink(old_zips)

  parts <- partition_family_files(all_files, max_bytes = MAX_PART)
  zip_names <- file.path(DIST_DIR, if (length(parts) == 1L)
    paste0(fid, "_models.zip") else sprintf("%s_models.part%d.zip", fid, seq_along(parts)))

  # Resolve both absolute paths BEFORE setwd(): afterwards MODEL_DIR is a
  # relative path that would resolve against itself ("models/models").
  abs_model_dir <- normalizePath(MODEL_DIR)
  failed <- FALSE
  for (k in seq_along(parts)) {
    zip_name  <- zip_names[k]
    message("  Creating: ", zip_name)
    abs_files <- normalizePath(parts[[k]])
    abs_zip   <- normalizePath(file.path(getwd(), zip_name), mustWork = FALSE)
    rel_files <- sub(paste0("^", abs_model_dir, .Platform$file.sep), "", abs_files)

    old_wd <- setwd(MODEL_DIR)
    zip::zip(zipfile = abs_zip, files = rel_files)
    setwd(old_wd)

    size_mb <- round(file.info(zip_name)$size / 1e6, 1)
    message("  Done: ", zip_name, " (", size_mb, " MB, ", length(abs_files), " files)")
    if (size_mb >= 100) {
      message("  TOO LARGE for GitHub (>= 100 MB) - lower MAX_PART; archive deleted")
      unlink(zip_name); failed <- TRUE; next
    }

    # Independent check on the bytes we are about to publish. The write-side
    # guards and the allowlist are promises made by code; this opens the
    # archive and verifies the property directly. A failure deletes the
    # archive so a dirty zip cannot be uploaded by mistake.
    audit <- audit_release_archive(zip_name, quiet = TRUE)
    if (!audit$ok) {
      message("  AUDIT FAILED for ", zip_name, " - archive deleted:")
      for (f in audit$findings) message("    ! ", f)
      unlink(zip_name); failed <- TRUE
      next
    }
    message("  Audit: PASS - no per-sample OSSL data in ", basename(zip_name))
  }
  if (failed) {
    # A family is published whole or not at all
    unlink(zip_names)
    audit_failures <- c(audit_failures, fid)
  }
}

if (length(audit_failures) > 0L) {
  message("\n=== DO NOT UPLOAD ===")
  message("The release audit failed for: ", paste(audit_failures, collapse = ", "))
  message("Those archives were deleted. Fix the artefacts and re-run before ",
          "publishing anything to Zenodo.")
  quit(status = 1L)
}

message("\n=== Upload checklist ===")
message("Every archive below was audited: no per-sample OSSL data, no file ",
        "outside the publishing allowlist.")
message("1. Copy into a clone of github.com/HugoMachadoRodrigues/autoSpectra-models:")
for (z in list.files(DIST_DIR, pattern = "_models(\\.part[0-9]+)?\\.zip$",
                     full.names = TRUE))
  if (sub("_models(\\.part[0-9]+)?\\.zip$", "", basename(z)) %in% FAMILIES)
    message("     ", normalizePath(z))
message("2. shasum -a 256 *.zip > SHA256SUMS ; commit ; push")
message("3. Create a GitHub release there; Zenodo archives it automatically")
message("4. Verify with: download_ossl_models(\"all\", overwrite = TRUE)")
