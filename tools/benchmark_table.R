#!/usr/bin/env Rscript
# Regenerate the README benchmark tables from the trained models' artefacts.
#
# The README tables must never be hand-maintained: they were wrong several
# times before (quick-mode numbers, models trained on rows with no spectra, and
# a Safanelli comparison against the wrong table). Emit them from the artefacts.
#
# Usage:  Rscript tools/benchmark_table.R            # markdown to stdout
#         Rscript tools/benchmark_table.R --check    # exit 1 if anything is missing
#         Rscript tools/benchmark_table.R --final-dir <dir>
#           read the held-out, window and interval tables from the published
#           models in <dir> (default models). The 10-fold CV tables always
#           come from models/ (cv10 runs).
#         Rscript tools/benchmark_table.R --lodo-dir <dir>
#           read the leave-one-dataset-out runs from <dir>/<family>/lodo/<prop>.json
#           (default models), as written by train_ossl.R --lodo --out-dir <dir>.

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE))
    pkgload::load_all(".", quiet = TRUE) else library(autoSpectra)
})

ARGS  <- commandArgs(trailingOnly = TRUE)
CHECK <- "--check" %in% ARGS
FINAL_DIR <- if (length(i <- which(ARGS == "--final-dir"))) ARGS[i[1] + 1L] else "models"
LODO_DIR  <- if (length(i <- which(ARGS == "--lodo-dir")))  ARGS[i[1] + 1L] else "models"
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
  agg_folds(pf, log_scale)
}

# One json per property (train_ossl.R --lodo): <dir>/<family>/lodo/<prop>.json
read_lodo <- function(family, prop, log_scale = FALSE, dir = LODO_DIR) {
  fp <- file.path(dir, family, "lodo", paste0(prop, ".json"))
  if (!file.exists(fp)) return(NULL)
  pf <- tryCatch(jsonlite::read_json(fp)$per_fold, error = function(e) NULL)
  agg_folds(pf, log_scale)
}

agg_folds <- function(pf, log_scale = FALSE) {
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
       k = sum(ok), n = sum(n[ok]), w_max = max(w))
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
cat("The published models. Cubist and PCA-VAE: the same 80/20 split, the same SNV spectra and the same 120 principal components (OSSL-style PCA: fitted on every library spectrum outside the test set, bands standardised); only the regression model differs.\n")
cat("LUCAS revisits of the same point are kept on one side of every split.\n\n")
cat("| Sensor | Property | PCA-VAE R² / RMSE / CCC | Cubist R² / RMSE / CCC |\n")
cat("|---|---|---|---|\n")
for (dom in DOMS) {
  cub <- paste0("OSSL_", dom, "_Cubist")
  pcv <- paste0("OSSL_", dom, "_PCAVAE")
  for (p in PROPS) {
    c <- read_info(cub, p, FINAL_DIR)$metrics; q <- read_info(pcv, p, FINAL_DIR)$metrics
    need(!is.null(c), paste0(cub, "/", p)); need(!is.null(q), paste0(pcv, "/", p))
    # Bold the best R² as shown (3 decimals), both cells on a tie
    best <- function(x) if (isTRUE(round(num(x$R2), 3) >=
                                   max(round(num(q$R2), 3), round(num(c$R2), 3)))) "**" else ""
    cell <- function(x) sprintf("%s%s%s / %s / %s", best(x), fmt(num(x$R2)), best(x),
                                fmt(num(x$RMSE), 2), fmt(num(x$CCC)))
    cat(sprintf("| %s | `%s` | %s | %s |\n", dom, p, cell(q), cell(c)))
  }
}

# ---- 3. PCA-VAE by instrument window ----------------------------------------

cat("\n### PCA-VAE by instrument window (R²)\n\n")
cat("Each window family is the PCA-VAE pipeline trained on one instrument window, on every OSSL spectrum that covers it; `predict_soil_agnostic()` sends each spectrum to the widest window it covers. Each family is scored on its own held-out 20%, so narrower windows are tested on more (and more varied) spectra.\n\n")
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
      m <- read_info(f, p, dir = FINAL_DIR)$metrics
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
      l  <- read_info(fid, p, FINAL_DIR)$latent
      cf <- l$conformal
      cat(sprintf("| %s | `%s` | %s | %s | %s | %s | %s |\n", dom, p,
                  if (alg == "PCAVAE") "PCA-VAE" else "Cubist",
                  pct(num(cf$test_coverage_68)), pct(num(cf$test_coverage_90)),
                  pct(num(cf$test_coverage_95)), pct(num(l$test_frac_out_of_domain))))
    }
  }
}

# ---- 5. Leave-one-dataset-out ----------------------------------------------

LODO_FAMS <- list(
  VisNIR = c(Cubist = "OSSL_VisNIR_Cubist", `PCA-VAE` = "OSSL_VisNIR_PCAVAE",
             `800–2500 nm` = "OSSL_VisNIR_PCAVAE_800",
             `1000–2500 nm` = "OSSL_VisNIR_PCAVAE_1000",
             `1350–2500 nm` = "OSSL_VisNIR_PCAVAE_1350"),
  MIR    = c(Cubist = "OSSL_MIR_Cubist", `PCA-VAE` = "OSSL_MIR_PCAVAE",
             `700–3800 cm⁻¹` = "OSSL_MIR_PCAVAE_700_3800"))
lodo_props <- function(fam) {
  f <- list.files(file.path(LODO_DIR, fam, "lodo"), pattern = "[.]json$")
  sub("[.]json$", "", f)
}

cat("\n### Leave-one-dataset-out, every family\n\n")
cat("Each fold withholds one whole OSSL contributing dataset; only properties measured in at least two datasets can be scored.",
    "R² is weighted by fold size, on the scale of the 10-fold table (`oc` and `n.tot` log1p; pH and clay original units).",
    "*Largest fold*: share of all test samples in the biggest withheld dataset. Best family per property in bold.\n\n")
for (dom in DOMS) {
  fams <- LODO_FAMS[[dom]]
  cat(sprintf("| %s | Folds | Largest fold | %s |\n", dom, paste(names(fams), collapse = " | ")))
  cat(sprintf("|---|---|---|%s\n", paste(rep("---|", length(fams)), collapse = "")))
  for (p in PROPS) {
    lg <- p %in% ossl_l1_log_transform
    rs <- lapply(fams, read_lodo, prop = p, log_scale = lg)
    for (k in names(fams)) need(!is.null(rs[[k]]), paste0(fams[[k]], " lodo/", p))
    r2 <- vapply(rs, function(r) if (is.null(r)) NA_real_ else r$r2_w, numeric(1))
    top <- if (any(is.finite(r2))) fmt(max(r2, na.rm = TRUE)) else ""
    cells <- vapply(r2, function(x) if (fmt(x) == top) paste0("**", fmt(x), "**") else fmt(x), "")
    r1 <- Filter(Negate(is.null), rs)[1][[1]]
    cat(sprintf("| `%s` | %s | %s | %s |\n", p, if (is.null(r1)) "—" else r1$k,
                if (is.null(r1)) "—" else pct(r1$w_max), paste(cells, collapse = " | ")))
  }
  cat("\n")
}

cat("All properties (on the scale of Safanelli's S2 Table: log1p except clay, silt, sand and pH):\n\n")
cat("| Sensor | Family | Properties | Median R² | R² > 0 | Best for | Median of median-fold R² | Best for (median fold) |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (dom in DOMS) {
  fams <- LODO_FAMS[[dom]]
  props <- Reduce(intersect, lapply(fams, lodo_props))
  need(length(props) > 0, paste0(dom, " lodo"))
  if (!length(props)) next
  get <- function(stat) sapply(fams, function(f) vapply(props, function(p) {
    r <- read_lodo(f, p, log_scale = p %in% ossl_l1_log_transform)
    if (is.null(r)) NA_real_ else r[[stat]] }, numeric(1)))
  W <- matrix(get("r2_w"), ncol = length(fams)); M <- matrix(get("r2_med"), ncol = length(fams))
  best <- function(X) tabulate(apply(X, 1, which.max), nbins = ncol(X))
  bw <- best(W); bm <- best(M)
  for (j in seq_along(fams))
    cat(sprintf("| %s | %s | %d | %s | %d | %d | %s | %d |\n", dom, names(fams)[j], length(props),
                fmt(stats::median(W[, j]), 2), sum(W[, j] > 0), bw[j],
                fmt(stats::median(M[, j]), 2), bm[j]))
}

if (CHECK) {
  if (length(missing)) {
    message("\nINCOMPLETE — do not paste this into the README. Missing: ",
            paste(unique(missing), collapse = ", "))
    quit(status = 1L)
  }
  message("\nComplete: every table is backed by an artefact.")
}
