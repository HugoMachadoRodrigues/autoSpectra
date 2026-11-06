# R/sensor_model.R - A PCA-VAE family tuned to one instrument
#
# The OSSL families are fitted on benchtop spectra. A different instrument
# (a handheld NeoSpectra, a field unit, a lab's own spectrometer) sees the same
# soils differently and often covers only part of the range, so it needs a
# model calibrated on its own labelled spectra. Which model is best depends on
# how many soils it has. On the OSSL NeoSpectra library (organic carbon, the
# same held-out US and African soils):
#
#   training soils   Cubist + PCA      PCA-VAE ensemble
#   150              0.66 / 0.71       0.60 / 0.25      (R2 USA / Africa)
#   400              0.73 / 0.81       0.71 / -1.10
#   1613             0.78 / 0.78       0.82 / 0.77
#
# With few soils the network can fail badly on soils unlike its training set;
# with many it is the more accurate. model = "auto" therefore uses Cubist below
# 1000 soils and the PCA-VAE from there on. Pretraining the PCA-VAE on OSSL cut
# to the instrument's range, then fine-tuning it, did not beat training on the
# instrument alone (1613 soils: 0.78 / 0.78), so it is optional.

#' Read a family description saved next to its models
#'
#' Families made by \code{train_sensor_model()} are not in
#' \code{model_registry}; their description is stored in
#' \code{<model_dir>/<family_id>/family.rds}.
#' @keywords internal
.custom_family <- function(id, model_dir) {
  fp <- file.path(model_dir, id, "family.rds")
  if (file.exists(fp)) readRDS(fp) else NULL
}

# Continue training each member on new data, keeping its stored predictor and
# target scaling (vae_fit() would recompute both and break the pretrained
# weights). Early stopping on the validation prediction loss keeps the best
# epoch; with no improvement the pretrained weights are kept.
.finetune_ensemble <- function(models, P, y, P_val, y_val, lr = 1e-4,
                               epochs = 200L, batch_size = 64L, patience = 15L) {
  keras <- reticulate::import("keras", convert = FALSE)
  lapply(models, function(m) {
    ym <- reticulate::py_to_r(m$y_mean$numpy())
    ys <- reticulate::py_to_r(m$y_sd$numpy())
    m$compile(optimizer = keras$optimizers$Adam(learning_rate = lr))
    m$fit(P, matrix((y - ym) / ys), epochs = as.integer(epochs),
          batch_size = as.integer(batch_size), verbose = 0L,
          validation_data = reticulate::tuple(P_val, matrix((y_val - ym) / ys)),
          callbacks = list(keras$callbacks$EarlyStopping(
            monitor = "val_y_loss", mode = "min", patience = as.integer(patience),
            restore_best_weights = TRUE)))
    m
  })
}

#' Train a model family for one instrument from its own labelled spectra
#'
#' Builds a model family for the instrument that measured \code{spectra}, on
#' the OSSL grid within the instrument's range and with the OSSL preprocessing.
#' The model is Cubist + PCA (as \code{train_cubist_pca()}) or a PCA + soilVAE
#' ensemble (as \code{train_pca_vae()}). \code{model = "auto"} picks Cubist
#' for fewer than 1000 distinct soils (\code{groups}) and the PCA-VAE from
#' 1000 on: on the OSSL NeoSpectra library the network was less accurate than
#' Cubist with 150 or 400 soils, and could fail on soils from another region,
#' but was more accurate with 1600.
#'
#' With \code{pretrain = TRUE} (PCA-VAE only) the ensemble is first trained on
#' the OSSL library cut to the instrument's range together with the
#' instrument's training spectra, then fine-tuned on the instrument's spectra.
#' On NeoSpectra this did not improve on training with the instrument's spectra
#' alone, so it is off by default.
#'
#' A grouped fifth of the samples is held out as a test set and the metrics on
#' it are saved with the model. Rows sharing a value of \code{groups} (for
#' example repeated scans of one soil) stay on the same side of every split;
#' give it whenever a soil was scanned more than once.
#'
#' The family is saved under \code{out_dir/family_id} and used with
#' \code{predict_soil(df, family_id, model_dir = out_dir)}. Calling the function
#' again for another property adds that property to the family.
#'
#' @param spectra Data frame with a \code{Soil_ID} column and one column per
#'   wavelength (nm) or wavenumber (cm-1), as for \code{predict_soil()}
#' @param y Lab values of \code{prop}, original units, one per row
#' @param prop Property name from \code{ossl_l1_properties} (sets the target
#'   transform and the OSSL column used for pretraining)
#' @param family_id Name of the new family, e.g. \code{"MyNeoSpectra"}
#' @param sensor_type \code{"visnir"} or \code{"mir"}
#' @param groups Optional sample label per row (see Details)
#' @param model \code{"auto"} (default), \code{"cubist"} or \code{"pca_vae"}
#' @param pretrain PCA-VAE only: pretrain on OSSL before fine-tuning (default
#'   \code{FALSE}); needs the OSSL files (downloaded from the official OSSL
#'   storage if missing and \code{download_if_missing = TRUE})
#' @param out_dir Root directory for the family's files
#' @param n_components,n_seeds PCA size and ensemble size
#' @param cache_dir,download_if_missing Passed to \code{ossl_prepare()}
#' @param seed Seed of the train/test split
#' @return Invisibly, the held-out test metrics (original units)
#' @export
train_sensor_model <- function(spectra, y, prop, family_id,
                               sensor_type  = c("visnir", "mir"),
                               groups       = NULL,
                               model        = c("auto", "cubist", "pca_vae"),
                               pretrain     = FALSE,
                               out_dir      = "models",
                               n_components = 120L,
                               n_seeds      = 5L,
                               cache_dir    = ossl_cache_dir(),
                               download_if_missing = TRUE,
                               seed         = 42L) {
  sensor_type <- match.arg(sensor_type)
  model <- match.arg(model)
  if (family_id %in% names(model_registry))
    stop("family_id '", family_id, "' is a built-in family; choose another name.")
  if (!prop %in% ossl_l1_properties) stop("Unknown property: ", prop)
  base_fam <- get_family(if (sensor_type == "visnir") "OSSL_VisNIR_Cubist" else "OSSL_MIR_Cubist")

  # The instrument's grid: the OSSL grid within the range it covers.
  wl  <- get_wavelengths(spectra, id_col = "Soil_ID")
  rng <- range(wl$wl)
  grid <- base_fam$wavegrid[base_fam$wavegrid >= rng[1] & base_fam$wavegrid <= rng[2]]
  if (length(grid) < 50L) stop("The spectra cover too few bands of the OSSL grid.")
  # Preprocess on the instrument's grid; `full` marks rows covering all of it
  # (prepare_model_input() zero-fills uncovered bands, so test the mask).
  prep <- function(M, src_wl) {
    R  <- resample_to_grid(M, src_wl = src_wl, target_wl = grid)
    mk <- mask_from_coverage(R)
    list(X = prepare_model_input(R, base_fam$preprocess, mask = mk)$X,
         full = rowSums(mk) == ncol(mk))
  }
  if (length(y) != nrow(spectra)) stop("y must have one value per row of spectra.")
  px <- prep(as.matrix(spectra[, wl$cols, drop = FALSE]), wl$wl)
  X  <- px$X
  yt <- apply_target_transform(y, prop)
  ok <- is.finite(yt) & px$full
  if (sum(ok) < 50L) stop("Need at least 50 complete labelled spectra; have ", sum(ok), ".")
  X <- X[ok, , drop = FALSE]; yt <- yt[ok]
  g <- if (is.null(groups)) seq_len(nrow(X)) else as.character(groups)[ok]

  old <- .custom_family(family_id, out_dir)
  if (!is.null(old) && !identical(old$wavegrid, grid))
    stop("Family '", family_id, "' already exists with a different spectral range.")
  want <- c(auto = NA, cubist = "cubist_pca", pca_vae = "pca_vae")[[model]]
  algo <- if (!is.null(old)) old$algorithm
          else if (!is.na(want)) want
          else if (length(unique(g)) < 1000L) "cubist_pca" else "pca_vae"
  if (!is.null(old) && !is.na(want) && !identical(want, algo))
    stop("Family '", family_id, "' already uses ", old$algorithm, ".")
  fam <- list(id = family_id,
              label = sprintf("%s \u2014 %s (%g-%g)", family_id,
                              if (algo == "cubist_pca") "Cubist + PCA" else "PCA + soilVAE ensemble",
                              min(grid), max(grid)),
              sensor_type = sensor_type, algorithm = algo,
              wavegrid = grid, preprocess = base_fam$preprocess,
              properties = union(old$properties, prop), source = "custom",
              pretrained_on_ossl = algo == "pca_vae" && isTRUE(pretrain))
  ncomp <- min(n_components, ncol(X), nrow(X) - 1L)
  # drop any cached copy of an older model with the same name
  rm(list = grep(paste0("::", family_id, "/", prop), ls(.autoSpectra_cache), value = TRUE,
                 fixed = TRUE), envir = .autoSpectra_cache)

  if (algo == "cubist_pca") {
    message("  [", family_id, "/", prop, "] Cubist + PCA (", length(unique(g)), " soils)")
    dir_create(file.path(out_dir, family_id))
    train_cubist_pca(X, yt, family_id, prop, out_dir = out_dir,
                     n_components = ncomp, groups = g)
    saveRDS(fam, file.path(out_dir, family_id, "family.rds"))
    info <- jsonlite::read_json(file.path(out_dir, family_id, "models",
                                          paste0(prop, "_metrics.json")))
    return(invisible(lapply(info$metrics, unlist)))
  }
  .need_soilVAE()

  sp <- .grouped_partition(nrow(X), c(0.8, 0.2), g, seed = seed)
  tr <- sp[[1]]; te <- sp[[2]]
  vs <- .grouped_partition(length(tr), c(0.85, 0.15), g[tr], seed = seed + 1L)
  fit <- tr[vs[[1]]]; val <- tr[vs[[2]]]

  if (isTRUE(pretrain)) {
    message("  Loading OSSL (", sensor_type, ") cut to ", rng[1], "-", rng[2], " ...")
    jd <- ossl_prepare(sensor_type = sensor_type, cache_dir = cache_dir,
                       properties = prop, download_if_missing = download_if_missing)
    owl <- get_wavelengths(jd, id_col = "Soil_ID")
    po  <- prep(as.matrix(jd[, owl$cols, drop = FALSE]), owl$wl)
    yo  <- apply_target_transform(jd[[prop]], prop)
    ko  <- is.finite(yo) & po$full
    Xo  <- po$X[ko, , drop = FALSE]; yo <- yo[ko]; rm(jd, po)
    pca <- .fit_pca(rbind(Xo, X[tr, , drop = FALSE]), n_components = ncomp)
    P   <- function(M) .project_pca(pca, M)
    message("  Pretraining on ", nrow(Xo), " OSSL + ", length(fit), " instrument spectra ...")
    models <- .fit_vae_ensemble(rbind(P(Xo), P(X[fit, , drop = FALSE])), c(yo, yt[fit]),
                                P(X[val, , drop = FALSE]), yt[val], n_seeds = n_seeds)
    rm(Xo)
    message("  Fine-tuning on the instrument's spectra ...")
    models <- .finetune_ensemble(models, P(X[fit, , drop = FALSE]), yt[fit],
                                 P(X[val, , drop = FALSE]), yt[val])
  } else {
    pca <- .fit_pca(X[tr, , drop = FALSE], n_components = ncomp)
    P   <- function(M) .project_pca(pca, M)
    models <- .fit_vae_ensemble(P(X[fit, , drop = FALSE]), yt[fit],
                                P(X[val, , drop = FALSE]), yt[val], n_seeds = n_seeds,
                                batch_size = 64L, epochs = 300L, patience = 20L)
  }

  # Held-out metrics, prediction intervals and applicability limit, as in
  # train_pca_vae(); the test set is split in halves (calibration / coverage).
  E    <- .ensemble_predict(models, P(X[te, , drop = FALSE]))
  yhat <- rowMeans(E)
  bounds <- target_bounds(yt[tr])
  mets <- metrics_original_scale(yt[te], yhat, prop, bounds = bounds)
  half <- .grouped_partition(length(te), c(0.5, 0.5), g[te], seed = seed + 2L)
  sdv  <- apply(E, 1, stats::sd); s0 <- stats::median(sdv[half[[1]]])
  if (!is.finite(s0) || s0 <= 0) s0 <- 1e-6
  score <- abs(yt[te] - yhat) / (sdv + s0)
  conf <- list(levels = .PCAVAE_LEVELS, s0 = s0, confidence = 0.68,
               alpha_grid = as.numeric(stats::quantile(score[half[[1]]], .PCAVAE_LEVELS,
                                                       names = FALSE)))
  conf$alpha_q <- conf$alpha_grid[which.min(abs(.PCAVAE_LEVELS - 0.68))]
  coverage <- vapply(c(0.68, 0.90, 0.95), function(lv)
    mean(score[half[[2]]] <= .conformal_alpha(conf, lv)), numeric(1))
  X_tr <- X[tr, , drop = FALSE]
  ucl_q <- as.numeric(stats::quantile(
    q_statistic(X_tr - .reconstruct_pca(pca, .project_pca(pca, X_tr))), 0.99, names = FALSE))

  # Save the family description and the models
  fam_dir <- file.path(out_dir, family_id)
  paths <- .pcavae_paths(out_dir, family_id, prop, n_seeds)
  dir_create(dirname(paths$base))
  saveRDS(fam, file.path(fam_dir, "family.rds"))
  for (i in seq_along(models)) .sv("vae_save")(models[[i]], paths$nets[i])
  saveRDS(list(pca = .strip_pca_scores(pca), n_seeds = n_seeds,
               soilVAE_version = as.character(utils::packageVersion("soilVAE"))),
          paths$bundle)
  save_scaler(0, 1, path_base = paths$base, bounds = bounds)
  save_metrics_json(paths$base, metrics = mets, conf_q = list(), latent = list(
    method = "q_statistic", ucl_q = ucl_q, ucl_q_method = "empirical_q99_train",
    n_components = ncomp, n_seeds = n_seeds, pretrained_on_ossl = isTRUE(pretrain),
    n_train = length(tr), n_test = length(te),
    conformal = c(conf, list(test_coverage_68 = coverage[1], test_coverage_90 = coverage[2],
                             test_coverage_95 = coverage[3],
                             n_coverage_eval = length(half[[2]])))))
  message(sprintf("  [%s/%s] held-out R2=%.3f RMSE=%.3f CCC=%.3f (n test=%d)",
                  family_id, prop, mets$R2, mets$RMSE, mets$CCC, length(te)))
  invisible(mets)
}
