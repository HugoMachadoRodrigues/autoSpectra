# tools/t2_limits.R - Add the Hotelling T2 limit to models trained before it
#
# Models trained with autoSpectra 0.5.0 record the Q-statistic limit but not
# the Hotelling T2 one. This script rebuilds each model's training rows from
# the official OSSL files (same inputs, same grouped 80/20 split, seed 42),
# checks their number against the n_train stored in the model, and writes the
# 99th percentile of their T2 into the model's metrics file as `ucl_t2`.
# Only that number is written; nothing per sample.
#
# Usage (from the package root):
#   Rscript tools/t2_limits.R --model-dir <dir> [--families A,B] [--dry-run]
#
# Families default to every Cubist / PCA-VAE family found in --model-dir.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

args <- commandArgs(TRUE)
opt  <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1L]
}
model_dir <- opt("--model-dir") %||% stop("--model-dir is required")
dry_run   <- "--dry-run" %in% args
fams <- opt("--families")
fams <- if (is.null(fams)) {
  list.dirs(model_dir, recursive = FALSE, full.names = FALSE)
} else {
  strsplit(fams, ",", fixed = TRUE)[[1]]
}

fam_info <- Filter(Negate(is.null), lapply(fams, function(fid) {
  f <- tryCatch(get_family(fid, model_dir), error = function(e) NULL)
  if (is.null(f) || !f$algorithm %in% c("cubist_pca", "pca_vae")) return(NULL)
  list(id = fid, fam = f)
}))
if (!length(fam_info)) stop("No Cubist / PCA-VAE family in ", model_dir)
# Training data per sensor and subset: the KSSL families (subset "kssl") were
# trained on the KSSL rows only, as train_ossl.R --subset kssl does.
data_key <- vapply(fam_info, function(x)
  paste(x$fam$sensor_type, x$fam$subset %||% "all"), character(1))

n_done <- 0L
for (key in unique(data_key)) {
  fi1 <- fam_info[[match(key, data_key)]]$fam
  jd  <- ossl_prepare(fi1$sensor_type, properties = ossl_l1_properties)
  if (identical(fi1$subset, "kssl"))
    jd <- jd[jd[["dataset.code_ascii_txt"]] %in% "KSSL.SSL", , drop = FALSE]
  for (fi in fam_info[data_key == key]) {
    pp  <- .family_training_inputs(fi$fam, jd)
    suf <- if (identical(fi$fam$algorithm, "pca_vae")) "_pcavae.rds" else "_cubist.rds"
    for (f in Sys.glob(file.path(model_dir, fi$id, "models", paste0("*", suf)))) {
      prop <- sub(suf, "", basename(f), fixed = TRUE)
      mpath <- file.path(dirname(f), paste0(prop, "_metrics.json"))
      if (!file.exists(mpath) || is.null(jd[[prop]])) next
      pca  <- readRDS(f)$pca
      sdev <- .pca_sdev(pca)
      if (is.null(sdev)) { message("skip ", fi$id, "/", prop, ": no component SDs"); next }
      y    <- apply_target_transform(jd[[prop]], prop)
      keep <- is.finite(y) & pp$has_spectra
      sp   <- .grouped_partition(sum(keep), c(0.80, 0.20), pp$groups[keep], seed = 42L)
      if (!identical(length(sp[[1]]), as.integer(pca$n_train))) {
        message("skip ", fi$id, "/", prop, ": rebuilt ", length(sp[[1]]),
                " training rows, model has ", pca$n_train)
        next
      }
      X  <- pp$X[keep, , drop = FALSE]
      t2 <- hotelling_t2(.project_pca(pca, X[sp[[1]], , drop = FALSE]), sdev)
      ucl <- as.numeric(stats::quantile(t2, 0.99, names = FALSE))
      above <- mean(hotelling_t2(.project_pca(pca, X[sp[[2]], , drop = FALSE]), sdev) > ucl)
      message(sprintf("%s/%s  n_train=%d  ucl_t2=%.1f  held-out above=%.3f",
                      fi$id, prop, length(sp[[1]]), ucl, above))
      if (!dry_run) {
        info <- jsonlite::read_json(mpath)
        info$latent$ucl_t2        <- ucl
        info$latent$ucl_t2_method <- "empirical_q99_train"
        info$latent$test_frac_above_t2 <- above
        jsonlite::write_json(info, mpath, auto_unbox = TRUE, digits = NA, pretty = TRUE)
      }
      n_done <- n_done + 1L
    }
    rm(pp); invisible(gc(verbose = FALSE))
  }
  rm(jd); invisible(gc(verbose = FALSE))
}
message("Done: ", n_done, " models", if (dry_run) " (dry run, nothing written)")
