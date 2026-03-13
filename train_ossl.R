# train_ossl.R
# Quick-start script to train soilVAE models on OSSL v1.2 data.
#
# This script:
#   1. Downloads OSSL VisNIR and MIR data from Google Cloud Storage (v1.2)
#   2. Trains soilVAE models for all OSSL L1 soil properties
#   3. Saves models to models/OSSL_VisNIR/ and models/OSSL_MIR/
#
# Preprocessing (canonical two-step SG):
#   VisNIR: Absorbance -> SG smooth (11,2) -> SG 1st derivative (11,2,1)
#   MIR   : SG smooth (11,2) -> SG 1st derivative (11,2,1)
#
# Run from the autoSpectra project root:
#   Rscript train_ossl.R
# -------------------------------------------------------------------------

# Install / load package (source if not installed)
if (!requireNamespace("autoSpectra", quietly = TRUE)) {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
} else {
  library(autoSpectra)
}

# keras / tensorflow are needed for training
if (!requireNamespace("keras", quietly = TRUE))
  stop("Install keras and tensorflow:\n",
       "  install.packages('keras')\n",
       "  keras::install_keras()")
library(keras)

# ---- Configuration -------------------------------------------------------

OUT_DIR    <- "models"        # where trained models are saved
CACHE_DIR  <- ossl_cache_dir() # OSSL data cache (~/.local/share/R/autoSpectra/...)
EPOCHS     <- 50L
BATCH_SIZE <- 64L

# Subset of properties to train (NULL = all ossl_l1_properties)
PROPERTIES <- NULL

# ---- Train VisNIR agnostic model ----------------------------------------
message("\n========== OSSL VisNIR ==========")
train_ossl_models(
  family_id           = "OSSL_VisNIR",
  out_dir             = OUT_DIR,
  cache_dir           = CACHE_DIR,
  properties          = PROPERTIES,
  download_if_missing = TRUE,
  epochs              = EPOCHS,
  batch_size          = BATCH_SIZE
)

# ---- Train MIR agnostic model -------------------------------------------
message("\n========== OSSL MIR ==========")
train_ossl_models(
  family_id           = "OSSL_MIR",
  out_dir             = OUT_DIR,
  cache_dir           = CACHE_DIR,
  properties          = PROPERTIES,
  download_if_missing = TRUE,
  epochs              = EPOCHS,
  batch_size          = BATCH_SIZE
)

message("\nAll done! Launch the app with:")
message("  autoSpectra::run_autoSpectra(model_dir = '", OUT_DIR, "')")
