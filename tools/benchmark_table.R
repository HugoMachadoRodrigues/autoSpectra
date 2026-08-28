#!/usr/bin/env Rscript
# Regenerate the README benchmark table from the trained models' metrics.
#
# The README table must never be hand-maintained: it was wrong twice before,
# once because it held quick-mode numbers and once because those numbers came
# from models trained on rows with no spectra. Emit it from the artefacts.
#
# Usage:  Rscript tools/benchmark_table.R            # markdown to stdout
#         Rscript tools/benchmark_table.R --check    # exit 1 if any model missing

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE))
    pkgload::load_all(".", quiet = TRUE) else library(autoSpectra)
})

CHECK <- "--check" %in% commandArgs(trailingOnly = TRUE)

# Safanelli et al. (2025) PLoS ONE 20(1):e0296545, Tables 4 (VisNIR) and 5 (MIR):
# Cubist + OSSL, mean over 20 instruments.
SAF <- list(
  VisNIR = list(oc       = c(R2 = 0.78, RMSE =  5.27, RPIQ = 0.65, CCC = 0.82),
                n.tot    = c(R2 = 0.81, RMSE =    NA, RPIQ = 2.26, CCC =   NA),
                ph.h2o   = c(R2 = 0.50, RMSE =  1.22, RPIQ = 2.18, CCC = 0.65),
                clay.tot = c(R2 = 0.25, RMSE = 13.35, RPIQ = 1.35, CCC = 0.38)),
  MIR    = list(oc       = c(R2 = 0.96, RMSE =  1.16, RPIQ = 2.96, CCC = 0.95),
                n.tot    = c(R2 = 0.96, RMSE =    NA, RPIQ = 4.98, CCC =   NA),
                ph.h2o   = c(R2 = 0.75, RMSE =  0.72, RPIQ = 2.92, CCC = 0.84),
                clay.tot = c(R2 = 0.62, RMSE =  8.82, RPIQ = 1.81, CCC = 0.74)))

PROPS <- c("oc", "n.tot", "ph.h2o", "clay.tot")

read_metric <- function(family, prop, field) {
  fp <- file.path("models", family, "models", paste0(prop, "_metrics.json"))
  if (!file.exists(fp)) return(NA_real_)
  m <- tryCatch(jsonlite::read_json(fp)$metrics, error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  as.numeric(unlist(m[[field]]))
}

fmt <- function(x, d = 3) if (!is.finite(x)) "—" else formatC(x, format = "f", digits = d)

missing <- character(0)

for (dom in c("VisNIR", "MIR")) {
  vae <- paste0("OSSL_", dom)
  cub <- paste0("OSSL_", dom, "_Cubist")

  cat(sprintf("\n### %s\n\n", dom))
  cat("| Property | soilVAE R² / RMSE / CCC | Cubist R² / RMSE / CCC | Safanelli R² / RMSE / CCC |\n")
  cat("|---|---|---|---|\n")

  for (p in PROPS) {
    for (f in c(vae, cub))
      if (!is.finite(read_metric(f, p, "R2"))) missing <- c(missing, paste0(f, "/", p))
    s <- SAF[[dom]][[p]]
    cat(sprintf("| `%s` | %s / %s / %s | %s / %s / %s | %s / %s / %s |\n", p,
        fmt(read_metric(vae, p, "R2")), fmt(read_metric(vae, p, "RMSE"), 2), fmt(read_metric(vae, p, "CCC")),
        fmt(read_metric(cub, p, "R2")), fmt(read_metric(cub, p, "RMSE"), 2), fmt(read_metric(cub, p, "CCC")),
        fmt(s[["R2"]], 2), fmt(s[["RMSE"]], 2), fmt(s[["CCC"]], 2)))
  }
}

cat("\n")
cat("> Metrics are test-set values in ORIGINAL units (log1p targets are\n")
cat("> back-transformed with expm1 before scoring), from a single 80/20 split.\n")
cat("> Safanelli's figures are means over 20 instruments on his own evaluation\n")
cat("> subset, NOT a global OSSL test set. For MIR in particular the populations\n")
cat("> differ: his reported RMSE and R2 imply an evaluation-set sd of about 5.8\n")
cat("> for oc against 11.6 in full OSSL, so the MIR columns are not a like-for-\n")
cat("> like comparison. The VisNIR populations do line up (implied sd 11.24 vs\n")
cat("> our 11.6).\n")

if (length(missing)) {
  message("\nMissing models: ", paste(unique(missing), collapse = ", "))
  if (CHECK) quit(status = 1L)
} else if (CHECK) {
  message("\nAll 8 family/property metrics present.")
}
