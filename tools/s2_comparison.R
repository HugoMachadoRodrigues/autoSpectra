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
# Runs that stored pooled out-of-fold metrics also report cor2 (squared
# correlation of all out-of-fold predictions), the R2 of the OSSL code.
#
# Cubist and PCA-VAE share their folds (seed 42), so their ten fold R2 values
# are paired: p_wilcoxon is a two-sided exact Wilcoxon signed-rank test of
# them (0.002 is the smallest p attainable with ten folds).
#
# Usage: Rscript tools/s2_comparison.R [--csv out.csv] [--kssl-dir D] [--ossl-dir D]
#   (markdown to stdout; the dirs default to models_kssl and models)

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(f, d) if (f %in% args) args[which(args == f) + 1L] else d
out_csv  <- opt("--csv", NULL)
ossl_dir <- opt("--ossl-dir", "models")
kssl_dir <- opt("--kssl-dir", "models_kssl")

s2 <- utils::read.csv("tools/safanelli2025_s2.csv", stringsAsFactors = FALSE)
code2prop <- setNames(names(.ossl_l1_columns), unname(.ossl_l1_columns))
s2$prop <- unname(code2prop[s2$property])
types <- list(
  mir_ossl            = c(dir = ossl_dir, dom = "MIR"),
  mir_kssl            = c(dir = kssl_dir, dom = "MIR"),
  visnir_ossl         = c(dir = ossl_dir, dom = "VisNIR"),
  visnir_kssl         = c(dir = kssl_dir, dom = "VisNIR"),
  nir.neospectra_ossl = c(dir = ossl_dir, dom = "NeoSpectra"))

read_cv <- function(dir, family, prop) {
  fp <- file.path(dir, family, "cv10", paste0(prop, ".json"))
  js <- if (file.exists(fp)) jsonlite::read_json(fp)
  pf <- if (!is.null(js)) js$per_fold else {
    fs <- file.path(dir, family, "cv10_summary.json")       # older layout
    if (file.exists(fs)) jsonlite::read_json(fs)$cv10[[prop]]$per_fold }
  list(js = js, pf = pf)
}
fold_val <- function(f, k, log_scale) { s <- if (log_scale) f$log_scale else f
  if (is.null(s[[k]])) NA_real_ else as.numeric(unlist(s[[k]])) }
fold_r2 <- function(dir, family, prop, log_scale)
  vapply(read_cv(dir, family, prop)$pf, fold_val, numeric(1), "R2", log_scale)

ours <- function(dir, family, prop, log_scale) {
  cv <- read_cv(dir, family, prop); js <- cv$js; pf <- cv$pf
  if (!length(pf)) return(c(R2 = NA, RMSE = NA, CCC = NA, cor2 = NA))
  g <- function(f, k) fold_val(f, k, log_scale)
  n <- vapply(pf, function(f) as.numeric(unlist(f$n_test)), numeric(1))
  r2 <- vapply(pf, g, numeric(1), "R2"); rm <- vapply(pf, g, numeric(1), "RMSE")
  cc <- vapply(pf, g, numeric(1), "CCC"); w <- n / sum(n)
  po <- js$pooled; if (log_scale && !is.null(po)) po <- po$log_scale
  c(R2 = sum(w * r2), RMSE = sqrt(sum(w * rm^2)), CCC = sum(w * cc),
    cor2 = if (is.null(po$cor2)) NA_real_ else as.numeric(po$cor2))
}

rows <- lapply(seq_len(nrow(s2)), function(i) {
  r <- s2[i, ]; ty <- types[[r$model]]
  lg <- r$prop %in% ossl_l1_log_transform
  cu <- ours(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_Cubist"), r$prop, lg)
  pv <- ours(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_PCAVAE"), r$prop, lg)
  fc <- fold_r2(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_Cubist"), r$prop, lg)
  fv <- fold_r2(ty[["dir"]], paste0("OSSL_", ty[["dom"]], "_PCAVAE"), r$prop, lg)
  pw <- if (length(fc) >= 5 && length(fc) == length(fv) && any(fv != fc))
    suppressWarnings(stats::wilcox.test(fv, fc, paired = TRUE, exact = TRUE)$p.value) else NA_real_
  data.frame(type = r$model, property = r$prop, scale = if (lg) "log1p" else "original",
             paper_R2 = r$R2, paper_RMSE = r$RMSE, cubist_R2 = cu[["R2"]],
             cubist_RMSE = cu[["RMSE"]], pcavae_R2 = pv[["R2"]], pcavae_RMSE = pv[["RMSE"]],
             cubist_cor2 = cu[["cor2"]], pcavae_cor2 = pv[["cor2"]], p_wilcoxon = pw)
})
tab <- do.call(rbind, rows)
if (!is.null(out_csv)) utils::write.csv(tab, out_csv, row.names = FALSE)

f3 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "—")
cat("| Type | Property | Scale | Paper R² | Our Cubist R² | PCA-VAE R² | PCA-VAE − Cubist | p (Wilcoxon) | Cubist cor² | PCA-VAE cor² |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf("| %s | `%s` | %s | %.2f | %s | %s | %s | %s | %s | %s |\n",
  type, property, scale, paper_R2, f3(cubist_R2), f3(pcavae_R2), f3(pcavae_R2 - cubist_R2),
  f3(p_wilcoxon), f3(cubist_cor2), f3(pcavae_cor2))))
done <- is.finite(tab$cubist_R2) & is.finite(tab$pcavae_R2)
message(sprintf("\n%d of %d models done on both sides. Our Cubist vs paper: median |dR2| = %s. PCA-VAE beats our Cubist in %d of %d.",
  sum(done), nrow(tab), f3(stats::median(abs(tab$cubist_R2 - tab$paper_R2)[is.finite(tab$cubist_R2)])),
  sum(tab$pcavae_R2[done] > tab$cubist_R2[done] + 0.005), sum(done)))
d <- tab$pcavae_R2 - tab$cubist_R2; sig <- done & is.finite(tab$p_wilcoxon) & tab$p_wilcoxon < 0.05
message(sprintf("Paired Wilcoxon (p < 0.05): PCA-VAE better in %d, tie in %d, worse in %d of %d.",
  sum(sig & d > 0), sum(done & !sig), sum(sig & d < 0), sum(done)))
for (ty in unique(tab$type)) { k <- done & tab$type == ty
  message(sprintf("  %-20s n=%2d  better=%2d tie=%2d worse=%d  mean R2 paper %.3f, Cubist %.3f, PCA-VAE %.3f",
    ty, sum(k), sum(k & sig & d > 0), sum(k & !sig), sum(k & sig & d < 0),
    mean(tab$paper_R2[k]), mean(tab$cubist_R2[k]), mean(tab$pcavae_R2[k]))) }
