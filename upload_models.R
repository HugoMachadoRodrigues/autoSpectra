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

  # Collect all files to include
  all_files <- list.files(src_dir, recursive = TRUE, full.names = TRUE)
  # Exclude .log / .tmp intermediates
  all_files <- all_files[!grepl("\\.(log|tmp)$", all_files)]

  zip_name <- file.path(DIST_DIR, paste0(fid, "_models.zip"))
  message("  Creating: ", zip_name)

  old_wd   <- setwd(MODEL_DIR)
  rel_files <- sub(paste0("^", normalizePath(MODEL_DIR), .Platform$file.sep),
                   "", normalizePath(all_files))
  zip::zip(zipfile = normalizePath(file.path(old_wd, zip_name), mustWork = FALSE),
           files   = rel_files)
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
