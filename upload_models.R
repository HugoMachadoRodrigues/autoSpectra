# upload_models.R — Package trained models for Zenodo upload
#
# Run AFTER train_ossl.R completes. Creates per-family zip archives for
# the autoSpectra Zenodo record (https://zenodo.org/record/19004686).
#
# Output:
#   dist/OSSL_VisNIR_models.zip            (masked soilVAE .keras + scalers)
#   dist/OSSL_MIR_models.zip               (masked soilVAE .keras + scalers)
#   dist/OSSL_VisNIR_Cubist_models.zip     (Cubist .rds + PCA loadings)
#   dist/OSSL_MIR_Cubist_models.zip        (Cubist .rds + PCA loadings)
#
# Each zip contains:
#   <family_id>/models/<prop>.keras / <prop>_cubist.rds / <prop>_cubist_error.rds
#   <family_id>/models/<prop>_scaler.rds / <prop>_metrics.json
#   <family_id>/metrics_summary.json
#
# Upload procedure:
#   1. Go to https://zenodo.org/record/19004686 → Edit → Upload files
#   2. Upload the zip archives
#   3. Save + Publish (new version of the same deposit — record ID stays the same)
#   4. Update .MODELS_ZENODO_RECORD in R/download_models.R if the record ID
#      ever changes (e.g. new deposit).
#
# Usage:
#   Rscript upload_models.R                      # all 4 families
#   Rscript upload_models.R OSSL_VisNIR          # one family
#   Rscript upload_models.R OSSL_VisNIR OSSL_MIR # two families

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

args      <- commandArgs(trailingOnly = TRUE)
ALL_FAM   <- c("OSSL_VisNIR", "OSSL_MIR", "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist")
FAMILIES  <- if (length(args) > 0) args else ALL_FAM
MODEL_DIR <- "models"
DIST_DIR  <- "dist"

dir_create(DIST_DIR)

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

  zip_name <- file.path(DIST_DIR, paste0(fid, "_models.zip"))
  message("  Creating: ", zip_name)

  # Resolve both absolute paths BEFORE setwd(): afterwards MODEL_DIR is a
  # relative path that would resolve against itself ("models/models").
  abs_model_dir <- normalizePath(MODEL_DIR)
  abs_files     <- normalizePath(all_files)
  abs_zip       <- normalizePath(file.path(getwd(), zip_name), mustWork = FALSE)

  rel_files <- sub(paste0("^", abs_model_dir, .Platform$file.sep), "", abs_files)

  old_wd <- setwd(MODEL_DIR)
  on.exit(setwd(old_wd), add = TRUE)
  zip::zip(zipfile = abs_zip, files = rel_files)
  setwd(old_wd)

  size_mb <- round(file.info(zip_name)$size / 1e6, 1)
  message("  Done: ", zip_name, " (", size_mb, " MB, ", length(all_files), " files)")
}

message("\n=== Upload checklist ===")
message("1. Go to: https://zenodo.org/record/19004686")
message("2. Click Edit -> New Version -> Upload files")
for (fid in FAMILIES) {
  zip_name <- file.path(DIST_DIR, paste0(fid, "_models.zip"))
  if (file.exists(zip_name))
    message("3. Upload: ", normalizePath(zip_name))
}
message("4. Save + Publish")
message("5. Verify with: download_ossl_models('OSSL_VisNIR', overwrite = TRUE)")
