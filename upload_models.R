# upload_models.R — Package trained models for Zenodo upload
#
# Run this AFTER train_ossl.R has completed successfully.
# Creates two zip archives ready for upload to the autoSpectra Zenodo record:
#
#   dist/OSSL_VisNIR_models.zip
#   dist/OSSL_MIR_models.zip
#
# Upload procedure:
#   1. Go to https://zenodo.org/record/19004686 → Edit → Upload files
#   2. Upload both zips
#   3. Save + Publish (creates a new version of the record)
#   4. Update .MODELS_ZENODO_RECORD in R/download_models.R if the
#      record ID changes (it stays the same for new versions of the same deposit)
#
# Usage:
#   Rscript upload_models.R
#   Rscript upload_models.R OSSL_VisNIR   # VisNIR only

library(autoSpectra)

args      <- commandArgs(trailingOnly = TRUE)
FAMILIES  <- if (length(args) > 0) args else c("OSSL_VisNIR", "OSSL_MIR")
MODEL_DIR <- getOption("autoSpectra.model_dir", "models")
DIST_DIR  <- "dist"

dir.create(DIST_DIR, showWarnings = FALSE)

for (fid in FAMILIES) {
  src_dir <- file.path(MODEL_DIR, fid)

  if (!dir.exists(src_dir)) {
    message("SKIP [", fid, "]: '", src_dir, "' not found.")
    message("  Run train_ossl.R first.")
    next
  }

  # Count trained models
  h5_files <- list.files(file.path(src_dir, "models"),
                          pattern = "\\.h5$", full.names = FALSE)
  message("[", fid, "]: ", length(h5_files), " trained models found.")

  if (length(h5_files) == 0) {
    message("  No .h5 files — skipping.")
    next
  }

  # List all files to include
  all_files <- list.files(src_dir, recursive = TRUE, full.names = TRUE)

  # Exclude large intermediate files if any
  all_files <- all_files[!grepl("\\.(log|tmp)$", all_files)]

  zip_name <- file.path(DIST_DIR, paste0(fid, "_models.zip"))

  message("  Creating: ", zip_name, " ...")
  old_wd <- setwd(MODEL_DIR)
  rel_files <- sub(paste0("^", normalizePath(MODEL_DIR), .Platform$file.sep),
                   "", normalizePath(all_files))
  zip::zip(zipfile = normalizePath(file.path(old_wd, zip_name),
                                    mustWork = FALSE),
            files   = rel_files)
  setwd(old_wd)

  size_mb <- round(file.info(zip_name)$size / 1e6, 1)
  message("  Done: ", zip_name, " (", size_mb, " MB)")
  message("  Files in archive: ", length(all_files))
}

message("\n=== Upload checklist ===")
message("1. Go to: https://zenodo.org/record/19004686")
message("2. Click 'Edit' -> 'Upload files'")
for (fid in FAMILIES) {
  zip_name <- file.path(DIST_DIR, paste0(fid, "_models.zip"))
  if (file.exists(zip_name))
    message("3. Upload: ", normalizePath(zip_name))
}
message("4. Save + Publish new version")
message("5. Verify download_ossl_models() works:")
message("   library(autoSpectra)")
message("   download_ossl_models('OSSL_VisNIR', overwrite = TRUE)")
