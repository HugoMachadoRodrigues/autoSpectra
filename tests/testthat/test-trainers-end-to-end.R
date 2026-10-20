# Run every trainer end to end on tiny data.
#
# The trainers are expensive on real data, so nothing exercised them and their
# failures surfaced only after hours of fitting. Three such defects reached the
# training path in a row, each one landing AFTER the fit:
#
#   * keras3::save_model() without overwrite = TRUE discarded every model on a
#     re-run ("File already exists").
#   * target_bounds(y_tr) referenced a variable that does not exist in
#     train_maskedVAE / train_soilVAE (the training targets there are y_tr_z,
#     z-scored), so both trainers would fit for their whole epoch budget and
#     then abort while saving the scaler. R CMD check found it as "no visible
#     binding for global variable 'y_tr'"; no test did.
#
# A few hundred samples and two epochs take seconds and walk the same code as a
# full run: fit, back-transform, score, save, reload.

tiny_spectra <- function(n = 240L, d = 40L, seed = 1L) {
  set.seed(seed)
  base <- sin(seq_len(d) / 6)
  X <- t(vapply(seq_len(n), function(i)
    0.4 + 0.15 * base * (1 + i / n) + stats::rnorm(d, sd = 0.01), numeric(d)))
  # a skewed, non-negative target driven by the spectra, so log1p applies
  y_raw <- exp(2 * rowMeans(X[, 5:15])) - 1 + abs(stats::rnorm(n, sd = 0.05))
  list(X = X, y_raw = y_raw, d = d)
}

skip_if_no_keras <- function() {
  skip_if_not_installed("keras3")
  skip_if_not(
    tryCatch({ autoSpectra:::.init_vae_python(); TRUE }, error = function(e) FALSE),
    "keras3 Python backend unavailable"
  )
}

expect_trained_artefacts <- function(out_dir, family, prop, model_file) {
  base <- file.path(out_dir, family, "models", prop)
  expect_true(file.exists(paste0(base, model_file)), info = model_file)
  expect_true(file.exists(paste0(base, "_metrics.json")))
  expect_true(file.exists(paste0(base, "_scaler.rds")))

  sc <- readRDS(paste0(base, "_scaler.rds"))
  # the calibration range the bounded back-transform depends on
  expect_true(is.finite(sc$bound_lo) && is.finite(sc$bound_hi))
  expect_lt(sc$bound_lo, sc$bound_hi)

  m <- jsonlite::read_json(paste0(base, "_metrics.json"))$metrics
  expect_true(is.finite(unlist(m$RMSE)))
  expect_true(is.finite(unlist(m$R2)))
}

test_that("train_maskedVAE fits, scores and saves end to end", {
  skip_if_no_keras()
  dat <- tiny_spectra()
  out <- file.path(tempdir(), "e2e_masked"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_no_error(
    train_maskedVAE(dat$X, apply_target_transform(dat$y_raw, "oc"),
                    family_id = "FAM", prop = "oc",
                    wavegrid = seq(400, by = 10, length.out = dat$d),
                    sensor_type = "visnir", out_dir = out,
                    epochs = 2L, batch_size = 32L, latent_dim = 4L,
                    kl_warmup_epochs = 1L, n_aug_copies = 1L, min_n = 20L)
  )
  expect_trained_artefacts(out, "FAM", "oc", ".keras")
})

test_that("train_soilVAE fits, scores and saves end to end", {
  skip_if_no_keras()
  dat <- tiny_spectra(seed = 2L)
  out <- file.path(tempdir(), "e2e_plain"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_no_error(
    train_soilVAE(dat$X, apply_target_transform(dat$y_raw, "oc"),
                  family_id = "FAM", prop = "oc", out_dir = out,
                  epochs = 2L, batch_size = 32L, latent_dim = 4L,
                  kl_warmup_epochs = 1L, min_n = 20L)
  )
  expect_trained_artefacts(out, "FAM", "oc", ".keras")
})

test_that("a VAE trainer can be re-run over its own output", {
  # The overwrite regression: the second run must not abort on save.
  skip_if_no_keras()
  dat <- tiny_spectra(seed = 3L)
  out <- file.path(tempdir(), "e2e_rerun"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  args <- list(dat$X, apply_target_transform(dat$y_raw, "oc"),
               family_id = "FAM", prop = "oc", out_dir = out,
               epochs = 1L, batch_size = 32L, latent_dim = 4L,
               kl_warmup_epochs = 1L, min_n = 20L)
  do.call(train_soilVAE, args)
  expect_no_error(do.call(train_soilVAE, args))
})

test_that("train_cubist_pca fits, scores and saves end to end", {
  skip_if_not_installed("Cubist")
  dat <- tiny_spectra(seed = 4L)
  out <- file.path(tempdir(), "e2e_cubist"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_no_error(
    train_cubist_pca(dat$X, apply_target_transform(dat$y_raw, "oc"),
                     family_id = "FAM", prop = "oc", out_dir = out,
                     n_components = 8L, committees_grid = c(1L, 2L),
                     k_inner = 2L, inner_subset_n = 200L, min_n = 20L)
  )
  expect_trained_artefacts(out, "FAM", "oc", "_cubist.rds")

  # and the saved bundle must not carry the training set (redistribution guard)
  b <- readRDS(file.path(out, "FAM", "models", "oc_cubist.rds"))
  expect_identical(b$model$data, "")
  expect_null(b$pca$scores)
})
