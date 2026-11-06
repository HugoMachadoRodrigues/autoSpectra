#!/usr/bin/env Rscript
# Compare our 10-fold CV results with Safanelli et al. (2025), S2 Table.
#
# tools/safanelli2025_s2.csv holds the paper's 141 rows (property code as the
# OSSL column name, model type, RMSE, bias, R2, CCC, RPIQ), transcribed from
# the supplementary PDF (journal.pone.0296545.s002). As in the paper, metrics
# are on the log1p scale except for clay, silt, sand and both pH.
#
# Our side: models/<family>/cv10/<prop>.json (KSSL types: models_kssl/...),
# written by `train_ossl.R --cv10`. Fold metrics are weighted by fold size.
#
# Usage: Rscript tools/s2_comparison.R [--csv out.csv]   (markdown to stdout)

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out_csv <- if ("--csv" %in% args) args[which(args == "--csv") + 1L] else NULL

s2 <- utils::read.csv("tools/safanelli2025_s2.csv", stringsAsFactors = FALSE)
code2prop <- setNames(names(.ossl_l1_columns), unname(.ossl_l1_columns))
s2$prop <- unname(code2prop[s2$property])
types <- list(
  mir_ossl            = c(dir = "models",      dom = "MIR"),
  mir_kssl            = c(dir = "models_kssl", dom = "MIR"),
  visnir_ossl         = c(dir = "models",      dom = "VisNIR"),
  visnir_kssl         = c(dir = "models_kssl", dom = "VisNIR"),
  nir.neospectra_ossl = c(dir = "models",      dom = "NeoSpectra"))

ours <- function(dir, family, prop, log_scale) {
  fp <- file.path(dir, family, "cv10", paste0(prop, ".json"))
  pf <- if (file.exists(fp)) jsonlite::read_json(fp)$per_fold else {
    fs <- file.path(dir, family, "cv10_summary.json")       # older layout
    if (file.exists(fs)) jsonlite::read_json(fs)$cv10[[prop]]$per_fold }
  if (!length(pf)) return(c(R2 = NA, RMSE = NA, CCC = NA))
  g <- function(f, k) { s <- if (log_scale) f$log_scale else f
    if (is.null(s[[k]])) NA_real_ else as.numeric(unlist(s[[k]])) }
  n <- vapply(pf, function(f) as.numeric(unlist(f$n_test)), numeric(1))
  r2 <- vapply(pf, g, numeric(1), "R2"); rm <- vapply(pf, g, numeric(1), "RMSE")
  cc <- vapply(pf, g, numeric(1), "CCC"); w <- n / sum(n)
  c(R2 = sum(w * r2), RMSE = sqrt(sum(w * rm^2)), CCC = sum(w * cc))
}

rows <- lapply(seq_len(nrow(s2)), function(i) {
  r <- s2[i, ]; ty <- types[[r$model]]
  lg <- r$prop %in% ossl_l1_log_transform
  cu <- ours(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_Cubist"), r$prop, lg)
  pv <- ours(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_PCAVAE"), r$prop, lg)
  data.frame(type = r$model, property = r$prop, scale = if (lg) "log1p" else "original",
             paper_R2 = r$R2, paper_RMSE = r$RMSE, cubist_R2 = cu[["R2"]],
             cubist_RMSE = cu[["RMSE"]], pcavae_R2 = pv[["R2"]], pcavae_RMSE = pv[["RMSE"]])
})
tab <- do.call(rbind, rows)
if (!is.null(out_csv)) utils::write.csv(tab, out_csv, row.names = FALSE)

f3 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "—")
cat("| Type | Property | Scale | Paper R² | Our Cubist R² | PCA-VAE R² | PCA-VAE − Cubist |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf("| %s | `%s` | %s | %.2f | %s | %s | %s |\n",
  type, property, scale, paper_R2, f3(cubist_R2), f3(pcavae_R2), f3(pcavae_R2 - cubist_R2))))
done <- is.finite(tab$cubist_R2) & is.finite(tab$pcavae_R2)
message(sprintf("\n%d of %d models done on both sides. Our Cubist vs paper: median |dR2| = %s. PCA-VAE beats our Cubist in %d of %d.",
  sum(done), nrow(tab), f3(stats::median(abs(tab$cubist_R2 - tab$paper_R2)[is.finite(tab$cubist_R2)])),
  sum(tab$pcavae_R2[done] > tab$cubist_R2[done] + 0.005), sum(done)))
