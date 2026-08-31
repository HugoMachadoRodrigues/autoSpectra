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

# LODO aggregation. The per-fold mean is misleading here: leave-one-dataset-out
# produces very few, very unequal folds (VisNIR has 4 datasets carrying spectra,
# with test sets from 50 to 40,175 samples), so a tiny fold swings the mean as
# hard as one 800x its size. Weight by fold size and report the fold count so
# the reader can judge.
read_lodo <- function(family, prop) {
  fp <- file.path("models", family, "cv10_summary.json")
  if (!file.exists(fp)) return(NULL)
  j <- tryCatch(jsonlite::read_json(fp), error = function(e) NULL)
  pf <- j$cv10[[prop]]$per_fold
  if (is.null(pf) || !length(pf)) return(NULL)
  num <- function(f, k) as.numeric(unlist(f[[k]]))
  r2 <- vapply(pf, num, numeric(1), "R2")
  n  <- vapply(pf, num, numeric(1), "n_test")
  ok <- is.finite(r2) & is.finite(n)
  if (!any(ok)) return(NULL)
  list(r2_w = sum(r2[ok] * n[ok]) / sum(n[ok]),
       r2_med = stats::median(r2[ok]),
       k = sum(ok))
}

missing <- character(0)
lodo_missing <- character(0)

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

have_lodo <- FALSE
for (dom in c("VisNIR", "MIR")) {
  cub <- paste0("OSSL_", dom, "_Cubist")
  rows <- lapply(PROPS, function(p) read_lodo(cub, p))
  names(rows) <- PROPS
  if (all(vapply(rows, is.null, logical(1)))) {
    lodo_missing <- c(lodo_missing, cub)
    next
  }
  if (!have_lodo) {
    cat("\n## Leave-one-dataset-out\n\n")
    cat("Each fold withholds an entire OSSL contributing dataset, so the model is\n")
    cat("scored on soil populations and instruments it never saw. This is the\n")
    cat("comparable protocol: Safanelli's figures are also means over held-out\n")
    cat("instruments. R² is weighted by fold size — LODO gives few, very unequal\n")
    cat("folds, and an unweighted mean lets a 50-sample fold outvote a 40,000-sample one.\n\n")
    have_lodo <- TRUE
  }
  cat(sprintf("\n### %s — Cubist under LODO\n\n", dom))
  cat("| Property | R² (n-weighted) | R² (median fold) | Folds | Safanelli R² |\n")
  cat("|---|---|---|---|---|\n")
  for (p in PROPS) {
    r <- rows[[p]]
    if (is.null(r)) {
      cat(sprintf("| `%s` | — | — | — | %s |\n", p, fmt(SAF[[dom]][[p]][["R2"]], 2)))
      lodo_missing <- c(lodo_missing, paste0(cub, "/", p))
    } else {
      cat(sprintf("| `%s` | %s | %s | %d | %s |\n", p,
          fmt(r$r2_w), fmt(r$r2_med), r$k, fmt(SAF[[dom]][[p]][["R2"]], 2)))
    }
  }
}

cat("\n")
cat("> **These columns are not a like-for-like comparison, and our numbers are\n")
cat("> not evidence of a better model.** Ours come from a single random 80/20\n")
cat("> split of the whole OSSL; Safanelli's are means over 20 instruments, each\n")
cat("> held out in turn. A random split leaves samples from the same dataset and\n")
cat("> region on both sides of the train/test line, so it is a materially easier\n")
cat("> evaluation than holding out whole instruments. Treat the gap as unexplained\n")
cat("> until a leave-one-dataset-out run under `--cv10 --lodo` measures it.\n")
cat("> \n")
cat("> Metrics are test-set values in ORIGINAL units: log1p targets are\n")
cat("> back-transformed with expm1 before scoring, so they are in the same units\n")
cat("> as Safanelli's tables.\n")

incomplete <- FALSE
if (length(missing)) {
  message("\nMissing single-split models: ", paste(unique(missing), collapse = ", "))
  incomplete <- TRUE
}
if (length(lodo_missing)) {
  message("\nMissing LODO results: ", paste(unique(lodo_missing), collapse = ", "))
  message("Run: Rscript train_ossl.R --only-cubist --cv10 --lodo --props ",
          paste(PROPS, collapse = ","))
  incomplete <- TRUE
}
if (CHECK) {
  if (incomplete) {
    message("\nINCOMPLETE — do not paste this into the README.")
    quit(status = 1L)
  }
  message("\nComplete: single-split and LODO present for all families.")
}
