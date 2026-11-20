#!/usr/bin/env Rscript
# Regenerate the README benchmark tables from the trained models' artefacts.
#
# The README tables must never be hand-maintained: they were wrong several
# times before (quick-mode numbers, models trained on rows with no spectra, and
# a Safanelli comparison against the wrong table). Emit them from the artefacts.
#
# Usage:  Rscript tools/benchmark_table.R            # markdown to stdout
#         Rscript tools/benchmark_table.R --check    # exit 1 if anything is missing
#         Rscript tools/benchmark_table.R --window-dir <dir>
#           read the instrument-window table (full-coverage PCA-VAE and the
#           window families, trained together) from <dir> (default models)

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE))
    pkgload::load_all(".", quiet = TRUE) else library(autoSpectra)
})

ARGS  <- commandArgs(trailingOnly = TRUE)
CHECK <- "--check" %in% ARGS
WINDOW_DIR <- if (length(i <- which(ARGS == "--window-dir"))) ARGS[i[1] + 1L] else "models"
PROPS <- c("oc", "n.tot", "ph.h2o", "clay.tot")
DOMS  <- c("VisNIR", "MIR")

# Safanelli et al. (2025) PLOS ONE 20(1):e0296545, S2 Table: Cubist trained on
# the OSSL subset, 10-fold cross-validation with refitting. Properties in
# ossl_l1_log_transform (here oc and n.tot) are reported on the log1p scale;
# clay and pH in original units. This is the like-for-like reference for a
# random partition of OSSL. (Tables 4 and 5 of the paper are a different test:
# fixed models scored on external ring-trial soils, see the README.)
S2 <- list(
  VisNIR = list(oc       = c(R2 = 0.94, RMSE = 0.24, CCC = 0.97),
                n.tot    = c(R2 = 0.92, RMSE = 0.08, CCC = 0.96),
                ph.h2o   = c(R2 = 0.90, RMSE = 0.43, CCC = 0.95),
                clay.tot = c(R2 = 0.82, RMSE = 6.61, CCC = 0.89)),
  MIR    = list(oc       = c(R2 = 0.99, RMSE = 0.12, CCC = 0.99),
                n.tot    = c(R2 = 0.97, RMSE = 0.06, CCC = 0.98),
                ph.h2o   = c(R2 = 0.91, RMSE = 0.40, CCC = 0.95),
                clay.tot = c(R2 = 0.90, RMSE = 5.64, CCC = 0.94)))

fmt <- function(x, d = 3) if (length(x) == 0 || !is.finite(x)) "—" else formatC(x, format = "f", digits = d)
pct <- function(x) if (length(x) == 0 || !is.finite(x)) "—" else sprintf("%.1f%%", 100 * x)
num <- function(x) if (is.null(x)) NA_real_ else as.numeric(unlist(x))

read_info <- function(family, prop, dir = "models") {
  fp <- file.path(dir, family, "models", paste0(prop, "_metrics.json"))
  if (!file.exists(fp)) return(NULL)
  tryCatch(jsonlite::read_json(fp), error = function(e) NULL)
}

# Fold-size-weighted aggregation of per-fold metrics. LODO folds are few and
# very unequal (VisNIR test sets from 50 to 40,175 samples), so an unweighted
# mean lets a tiny fold outvote a huge one. RMSE combines in quadrature.
read_folds <- function(family, file, prop, log_scale = FALSE) {
  fp <- file.path("models", family, file)
  if (!file.exists(fp)) return(NULL)
  pf <- tryCatch(jsonlite::read_json(fp)$cv10[[prop]]$per_fold, error = function(e) NULL)
  if (is.null(pf) || !length(pf)) return(NULL)
  get <- function(f, k) {
    src <- if (log_scale) f$log_scale else f
    if (is.null(src)) NA_real_ else num(src[[k]])
  }
  r2   <- vapply(pf, get, numeric(1), "R2")
  rmse <- vapply(pf, get, numeric(1), "RMSE")
  ccc  <- vapply(pf, get, numeric(1), "CCC")
  n    <- vapply(pf, function(f) num(f$n_test), numeric(1))
  ok <- is.finite(r2) & is.finite(n)
  if (!any(ok)) return(NULL)
  w <- n[ok] / sum(n[ok])
  list(r2_w = sum(r2[ok] * w), r2_med = stats::median(r2[ok]),
       rmse_w = sqrt(sum(w * rmse[ok]^2)), ccc_w = sum(ccc[ok] * w),
       k = sum(ok), n = sum(n[ok]))
}

missing <- character(0)
need <- function(ok, what) if (!ok) missing <<- c(missing, what)

# ---- 1. Comparison with Safanelli (random 10-fold CV) ----------------------

cat("### Cubist vs Safanelli et al. (2025): 10-fold cross-validation\n\n")
cat("Same protocol on both sides: Cubist + PCA(120) + SNV, random 10-fold CV on OSSL v1.2.\n")
cat("`oc` and `n.tot` are on the log1p scale, as the paper reports them; pH and clay in original units.\n\n")
cat("| Sensor | Property | Scale | autoSpectra R² / RMSE / CCC | Safanelli S2 R² / RMSE / CCC |\n")
cat("|---|---|---|---|---|\n")
for (dom in DOMS) {
  cub <- paste0("OSSL_", dom, "_Cubist")
  for (p in PROPS) {
    lg <- p %in% ossl_l1_log_transform
    r  <- read_folds(cub, "cv10_summary.json", p, log_scale = lg)
    need(!is.null(r), paste0(cub, " cv10/", p))
    s  <- S2[[dom]][[p]]
    cat(sprintf("| %s | `%s` | %s | %s / %s / %s | %s / %s / %s |\n", dom, p,
                if (lg) "log1p" else "original",
                fmt(r$r2_w), fmt(r$rmse_w, 2), fmt(r$ccc_w),
                fmt(s[["R2"]], 2), fmt(s[["RMSE"]], 2), fmt(s[["CCC"]], 2)))
  }
}

# ---- 2. Single split: Cubist vs PCA-VAE ------------------------------------

cat("\n### Held-out test set (original units)\n\n")
cat("Cubist and PCA-VAE: the same 80/20 split, the same SNV spectra and the same 120 principal components; only the regression model differs.\n")
cat("LUCAS revisits of the same point are kept on one side of every split.\n\n")
cat("| Sensor | Property | PCA-VAE R² / RMSE / CCC | Cubist R² / RMSE / CCC |\n")
cat("|---|---|---|---|\n")
for (dom in DOMS) {
  cub <- paste0("OSSL_", dom, "_Cubist")
  pcv <- paste0("OSSL_", dom, "_PCAVAE")
  for (p in PROPS) {
    c <- read_info(cub, p)$metrics; q <- read_info(pcv, p)$metrics
    need(!is.null(c), paste0(cub, "/", p)); need(!is.null(q), paste0(pcv, "/", p))
    best <- function(x) if (isTRUE(num(x$R2) >= max(num(q$R2), num(c$R2)))) "**" else ""
    cell <- function(x) sprintf("%s%s%s / %s / %s", best(x), fmt(num(x$R2)), best(x),
                                fmt(num(x$RMSE), 2), fmt(num(x$CCC)))
    cat(sprintf("| %s | `%s` | %s | %s |\n", dom, p, cell(q), cell(c)))
  }
}

# ---- 3. PCA-VAE by instrument window ----------------------------------------

cat("\n### PCA-VAE by instrument window (R²)\n\n")
cat("Each window family is the PCA-VAE pipeline trained on one instrument window, on every OSSL spectrum that covers it; `predict_soil_agnostic()` sends each spectrum to the widest window it covers. Each family is scored on its own held-out 20%, so narrower windows are tested on more (and more varied) spectra. All the families in this table, the full-range one included, use the OSSL-style PCA (`train_ossl.R --pca-fit library --pca-scale`).\n\n")
for (dom in DOMS) {
  fams <- autoSpectra:::.AGNOSTIC_WINDOWS[[tolower(dom)]]
  hdr  <- vapply(fams, function(f) {
    r <- range(get_family(f)$wavegrid)
    sprintf("%g–%g%s", r[1], r[2], if (dom == "MIR") " cm⁻¹" else " nm")
  }, character(1))
  cat(sprintf("| %s | %s |\n", dom, paste(hdr, collapse = " | ")))
  cat(sprintf("|---|%s\n", paste(rep("---|", length(fams)), collapse = "")))
  for (p in PROPS) {
    cells <- vapply(fams, function(f) {
      m <- read_info(f, p, dir = WINDOW_DIR)$metrics
      need(!is.null(m), paste0(f, "/", p))
      fmt(num(m$R2), 3)
    }, character(1))
    cat(sprintf("| `%s` | %s |\n", p, paste(cells, collapse = " | ")))
  }
  cat("\n")
}

# ---- 4. Intervals and applicability domain ----------------------------------

cat("### Prediction intervals and applicability domain (held-out rows)\n\n")
cat("| Sensor | Property | Family | 68% interval covers | 90% covers | 95% covers | Flagged out of domain |\n")
cat("|---|---|---|---|---|---|---|\n")
for (dom in DOMS) {
  for (alg in c("PCAVAE", "Cubist")) {
    fid <- paste0("OSSL_", dom, "_", alg)
    for (p in PROPS) {
      l  <- read_info(fid, p)$latent
      cf <- l$conformal
      cat(sprintf("| %s | `%s` | %s | %s | %s | %s | %s |\n", dom, p,
                  if (alg == "PCAVAE") "PCA-VAE" else "Cubist",
                  pct(num(cf$test_coverage_68)), pct(num(cf$test_coverage_90)),
                  pct(num(cf$test_coverage_95)), pct(num(l$test_frac_out_of_domain))))
    }
  }
}

# ---- 5. Leave-one-dataset-out ----------------------------------------------

cat("\n### Cubist under leave-one-dataset-out (original units)\n\n")
cat("Each fold withholds one whole OSSL contributing dataset. R² and CCC are weighted by fold size, RMSE pooled.\n\n")
cat("| Sensor | Property | R² | RMSE | CCC | R² median fold | Folds |\n")
cat("|---|---|---|---|---|---|---|\n")
for (dom in DOMS) {
  cub <- paste0("OSSL_", dom, "_Cubist")
  for (p in PROPS) {
    r <- read_folds(cub, "lodo_summary.json", p)
    need(!is.null(r), paste0(cub, " lodo/", p))
    cat(sprintf("| %s | `%s` | %s | %s | %s | %s | %s |\n", dom, p,
                fmt(r$r2_w), fmt(r$rmse_w, 2), fmt(r$ccc_w), fmt(r$r2_med),
                if (is.null(r)) "—" else r$k))
  }
}

if (CHECK) {
  if (length(missing)) {
    message("\nINCOMPLETE — do not paste this into the README. Missing: ",
            paste(unique(missing), collapse = ", "))
    quit(status = 1L)
  }
  message("\nComplete: every table is backed by an artefact.")
}
