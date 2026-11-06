#!/usr/bin/env Rscript
# Run 10-fold CV for Cubist and PCA-VAE on every model of Safanelli et al.
# (2025, S2 Table), one model type at a time. Resumable: finished properties
# are skipped (models/<family>/cv10/<prop>.json).
#
# Usage: Rscript tools/run_paper_replication.R visnir_ossl mir_ossl nir.neospectra_ossl
#        Rscript tools/run_paper_replication.R mir_kssl visnir_kssl
# Logs go to <out-dir>/paper_replication_<type>_<family>.log

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Give one or more model types: visnir_ossl, mir_ossl, nir.neospectra_ossl, mir_kssl, visnir_kssl")
s2 <- utils::read.csv("tools/safanelli2025_s2.csv", stringsAsFactors = FALSE)
env <- new.env(); sys.source("R/registry.R", envir = env)
code2prop <- setNames(names(env$.ossl_l1_columns), unname(env$.ossl_l1_columns))
spec <- list(visnir_ossl = c("VisNIR", "ossl"), mir_ossl = c("MIR", "ossl"),
             nir.neospectra_ossl = c("NeoSpectra", "ossl"),
             mir_kssl = c("MIR", "kssl"), visnir_kssl = c("VisNIR", "kssl"))
for (ty in args) {
  sp <- spec[[ty]]; if (is.null(sp)) stop("Unknown type: ", ty)
  props <- unname(code2prop[s2$property[s2$model == ty]])
  out_dir <- if (sp[2] == "kssl") "models_kssl" else "models"
  dir.create(out_dir, showWarnings = FALSE)
  for (alg in c("Cubist", "PCAVAE")) {
    fam <- paste0("OSSL_", sp[1], "_", alg)
    log <- file.path(out_dir, sprintf("paper_replication_%s_%s.log", ty, alg))
    message(format(Sys.time()), "  ", ty, " / ", fam, " : ", length(props), " properties -> ", log)
    st <- system2("Rscript", c("train_ossl.R", fam, "--cv10", "--subset", sp[2], "--out-dir", out_dir,
                              "--props", paste(props, collapse = ",")), stdout = log, stderr = log)
    if (st != 0) message("  exit status ", st, " - see ", log)
  }
}
