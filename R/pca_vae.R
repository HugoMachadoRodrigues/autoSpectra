# R/pca_vae.R - PCA + soilVAE ensemble ("PCA-VAE") family
#
# Same input as the Cubist family (SNV spectra -> 120 principal components,
# fitted on the training rows) and the same 80/20 split, so the two families
# differ only in the regression model. The model is an ensemble of supervised
# VAEs from the soilVAE package (>= 0.2.0), one per seed; the prediction is
# their mean.
#
# Why PCs rather than the raw spectrum: with 1000+ strongly correlated bands a
# dense network has far more weights than the data can pin down; the 120 PCs
# keep ~99.9% of the spectral variance and remove most of the noise. Why an
# ensemble: single networks differ from seed to seed; their mean is more
# accurate and steadier than any one of them.

.PCAVAE_LEVELS <- seq(0.50, 0.99, by = 0.01)

.need_soilVAE <- function() {
  if (!requireNamespace("soilVAE", quietly = TRUE) ||
      utils::packageVersion("soilVAE") < "0.2.0")
    stop("The PCA-VAE family needs the 'soilVAE' package (>= 0.2.0): ",
         "install.packages('soilVAE').", call. = FALSE)
  if (!isTRUE(soilVAE::vae_available()))
    stop("soilVAE needs Python with TensorFlow and Keras: run soilVAE::vae_install().",
         call. = FALSE)
  invisible(TRUE)
}

.pcavae_paths <- function(model_dir, family_id, prop, n_seeds = NULL) {
  base <- file.path(model_dir, family_id, "models", prop)
  list(base = base, bundle = paste0(base, "_pcavae.rds"),
       nets = if (!is.null(n_seeds)) sprintf("%s_vae%d.keras", base, seq_len(n_seeds)))
}

# Fit one ensemble on PC scores. `P_fit`/`P_val` are PC scores, `y_*` the
# targets on the model scale. Returns the list of fitted soilVAE models.
.fit_vae_ensemble <- function(P_fit, y_fit, P_val, y_val, n_seeds = 5L,
                              hidden_enc = c(256L, 128L), hidden_dec = c(128L, 256L),
                              latent_dim = 32L, epochs = 200L, batch_size = 256L,
                              patience = 15L, verbose = FALSE) {
  .need_soilVAE()
  lapply(seq_len(n_seeds), function(s) {
    try(reticulate::import("keras")$backend$clear_session(), silent = TRUE)
    m <- soilVAE::vae_build(ncol(P_fit), hidden_enc = hidden_enc, hidden_dec = hidden_dec,
                            latent_dim = latent_dim, seed = s)
    # "band" scaling: each PC has its own variance, which carries no meaning
    # for the network, so every PC is standardised separately.
    soilVAE::vae_fit(m, P_fit, y_fit, P_val, y_val, epochs = epochs,
                     batch_size = batch_size, patience = patience, scale_x = "band")
    if (verbose) message(sprintf("    seed %d/%d done", s, n_seeds))
    m
  })
}

# Matrix of predictions, one column per ensemble member.
.ensemble_predict <- function(models, P) {
  vapply(models, function(m) soilVAE::vae_predict(m, P), numeric(nrow(P)))
}

#' Train a PCA-VAE ensemble for one soil property
#'
#' Uses the same 80/20 split (grouped, seed 42) and the same 120-component PCA
#' as \code{train_cubist_pca()}, then fits \code{n_seeds} supervised VAEs from
#' the \pkg{soilVAE} package on the PC scores. A grouped tenth of the training
#' rows is held out for early stopping. The held-out 20\% is split in two
#' halves: one calibrates the prediction intervals, the other measures their
#' coverage; point metrics use all of it.
#'
#' Prediction intervals are normalised conformal intervals whose width follows
#' the spread of the ensemble members, so uncertain samples get wider
#' intervals.
#'
#' @inheritParams train_cubist_pca
#' @param n_seeds Number of ensemble members (default 5)
#' @param hidden_enc,hidden_dec,latent_dim Network sizes passed to
#'   \code{soilVAE::vae_build()}
#' @param epochs,batch_size,patience Training settings passed to
#'   \code{soilVAE::vae_fit()}
#' @return Invisible \code{TRUE}, or \code{NULL} when skipped
#' @export
train_pca_vae <- function(X, y, family_id, prop,
                          out_dir      = "models",
                          n_components = 120L,
                          n_seeds      = 5L,
                          hidden_enc   = c(256L, 128L),
                          hidden_dec   = c(128L, 256L),
                          latent_dim   = 32L,
                          epochs       = 200L,
                          batch_size   = 256L,
                          patience     = 15L,
                          min_n        = 50L,
                          groups       = NULL) {
  .need_soilVAE()
  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }
  if (is.null(groups)) groups <- seq_len(n)

  # Same partitions as train_cubist_pca(): identical test rows and halves.
  sp   <- .grouped_partition(n, c(0.80, 0.20), groups, seed = 42L)
  X_tr <- X[sp[[1]], , drop = FALSE]; y_tr <- y[sp[[1]]]; g_tr <- groups[sp[[1]]]
  X_te <- X[sp[[2]], , drop = FALSE]; y_te <- y[sp[[2]]]
  half   <- .grouped_partition(length(y_te), c(0.5, 0.5), groups[sp[[2]]], seed = 44L)
  cal_ix <- half[[1]]; eval_ix <- half[[2]]

  pca  <- .fit_pca(X_tr, n_components = n_components)
  P_tr <- pca$scores
  P_te <- .project_pca(pca, X_te)

  vs <- .grouped_partition(length(y_tr), c(0.9, 0.1), g_tr, seed = 7L)
  models <- .fit_vae_ensemble(P_tr[vs[[1]], , drop = FALSE], y_tr[vs[[1]]],
                              P_tr[vs[[2]], , drop = FALSE], y_tr[vs[[2]]],
                              n_seeds = n_seeds, hidden_enc = hidden_enc,
                              hidden_dec = hidden_dec, latent_dim = latent_dim,
                              epochs = epochs, batch_size = batch_size,
                              patience = patience)

  E_te    <- .ensemble_predict(models, P_te)
  yhat_te <- rowMeans(E_te)
  bounds  <- target_bounds(y_tr)
  mets    <- metrics_original_scale(y_te, yhat_te, prop, bounds = bounds)

  # Normalised conformal intervals: score = |residual| / (ensemble SD + s0).
  # s0 (median ensemble SD on the calibration half) keeps the width away from
  # zero where the members happen to agree.
  sd_te <- apply(E_te, 1, stats::sd)
  s0    <- stats::median(sd_te[cal_ix])
  if (!is.finite(s0) || s0 <= 0) s0 <- 1e-6
  delta <- sd_te + s0
  score <- abs(y_te - yhat_te) / delta
  alpha_grid <- as.numeric(stats::quantile(score[cal_ix], .PCAVAE_LEVELS, names = FALSE))
  conf <- list(levels = .PCAVAE_LEVELS, alpha_grid = alpha_grid, s0 = s0,
               confidence = 0.68, alpha_q = alpha_grid[which.min(abs(.PCAVAE_LEVELS - 0.68))])
  coverage <- vapply(c(0.68, 0.90, 0.95), function(lv)
    mean(score[eval_ix] <= .conformal_alpha(conf, lv)), numeric(1))

  # Applicability domain: Q-statistic on the same PCA, as for Cubist.
  E_tr  <- X_tr - .reconstruct_pca(pca, pca$scores)
  ucl_q <- as.numeric(stats::quantile(q_statistic(E_tr), 0.99, names = FALSE))
  rm(E_tr)
  Q_te  <- predict_q_statistic(pca, X_te[eval_ix, , drop = FALSE], ucl_q)

  # Save. REDISTRIBUTION GUARD: the PCA keeps only its basis (no scores); the
  # networks hold weights and per-PC means/SDs, no samples.
  paths <- .pcavae_paths(out_dir, family_id, prop, n_seeds)
  dir_create(dirname(paths$base))
  for (i in seq_along(models)) soilVAE::vae_save(models[[i]], paths$nets[i])
  saveRDS(list(pca = .strip_pca_scores(pca), n_seeds = n_seeds,
               soilVAE_version = as.character(utils::packageVersion("soilVAE"))),
          paths$bundle)
  save_scaler(0, 1, path_base = paths$base, bounds = bounds)
  save_metrics_json(paths$base, metrics = mets, conf_q = list(), latent = list(
    method                  = "q_statistic",
    ucl_q                   = ucl_q,
    ucl_q_method            = "empirical_q99_train",
    test_frac_out_of_domain = mean(!Q_te$in_domain),
    n_components            = n_components,
    var_explained           = sum(pca$var_explained),
    n_seeds                 = n_seeds,
    n_train                 = length(y_tr),
    n_test                  = length(y_te),
    conformal = c(conf, list(test_coverage_68 = coverage[1],
                             test_coverage_90 = coverage[2],
                             test_coverage_95 = coverage[3],
                             n_coverage_eval  = length(eval_ix)))))

  message(sprintf("  [%s/%s] RMSE=%.3f Bias=%.3f R2=%.3f RPIQ=%.2f CCC=%.3f seeds=%d n=%d",
                  family_id, prop, mets$RMSE, mets$Bias, mets$R2, mets$RPIQ,
                  mets$CCC, n_seeds, n))
  invisible(TRUE)
}

#' Load a saved PCA-VAE ensemble
#'
#' @param family_id Family ID
#' @param prop Soil property name
#' @param model_dir Root model directory
#' @return List with \code{pca}, \code{models} (soilVAE models) and
#'   \code{conformal} (interval calibration, or \code{NULL})
#' @export
load_pca_vae <- function(family_id, prop, model_dir = "models") {
  .need_soilVAE()
  p <- .pcavae_paths(model_dir, family_id, prop)
  if (!file.exists(p$bundle)) stop("PCA-VAE bundle not found: ", p$bundle)
  b <- readRDS(p$bundle)
  nets <- .pcavae_paths(model_dir, family_id, prop, b$n_seeds)$nets
  if (!all(file.exists(nets))) stop("Missing ensemble members for ", family_id, "/", prop)
  mp <- paste0(p$base, "_metrics.json")
  conf <- if (file.exists(mp)) jsonlite::read_json(mp)$latent$conformal
  list(pca = b$pca, models = lapply(nets, soilVAE::vae_load), conformal = conf)
}

.predict_one_pca_vae <- function(family_id, prop, X_proc, model_dir,
                                 with_pi = FALSE, confidence = 0.68) {
  if (!file.exists(.pcavae_paths(model_dir, family_id, prop)$bundle))
    return(list(y = rep(NA_real_, nrow(X_proc))))
  key <- paste0(.cache_key(model_dir, family_id, prop), "::pcavae")
  if (!exists(key, envir = .autoSpectra_cache, inherits = FALSE))
    assign(key, load_pca_vae(family_id, prop, model_dir), envir = .autoSpectra_cache)
  b <- get(key, envir = .autoSpectra_cache, inherits = FALSE)

  E    <- .ensemble_predict(b$models, .project_pca(b$pca, X_proc))
  if (is.null(dim(E))) E <- matrix(E, nrow = 1L)
  yhat <- rowMeans(E)
  res  <- list(y = yhat)
  if (with_pi && !is.null(b$conformal)) {
    a     <- .conformal_alpha(b$conformal, confidence)
    delta <- apply(E, 1, stats::sd) + as.numeric(b$conformal$s0)
    res$pi_lo <- yhat - a * delta
    res$pi_hi <- yhat + a * delta
  }
  res
}
