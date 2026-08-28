# run.R — Launch full OSSL training for both sensor-agnostic families.
#
# OSSL_VisNIR: Masked VAE — one model for all VisNIR sensors (400-2500 nm)
# OSSL_MIR:   Standard VAE — one model for all MIR sensors (600-4000 cm-1)
#
# Usage:
#   Rscript run.R                      # both families, all 34 properties
#   Rscript run.R --props oc,clay.tot  # sentinel properties only
#   Rscript run.R --quick              # quick validation run

message("run.R: starting autoSpectra full training ...")
flush.console()

args_raw <- c(
  "OSSL_VisNIR", "OSSL_MIR",
  "--epochs", "80",
  "--seed",   "42",
  "--workers", as.character(max(1L, parallel::detectCores(logical = FALSE) - 2L))
)

# Inherit any extra flags from the command line (e.g. --quick --props oc)
extra <- commandArgs(trailingOnly = TRUE)
extra <- extra[!extra %in% c("OSSL_VisNIR", "OSSL_MIR")]  # avoid duplication
args_raw <- c(args_raw, extra)

message("run.R: args = ", paste(args_raw, collapse = " "))
flush.console()

source("train_ossl.R")
