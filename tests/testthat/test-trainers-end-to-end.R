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

  info <- jsonlite::read_json(paste0(base, "_metrics.json"))
  # metrics and aggregate statistics only: no per-band training ranges
  expect_null(info$feat_range)
  m <- info$metrics
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

skip_if_no_soilVAE <- function() {
  skip_if_not_installed("soilVAE", minimum_version = "0.2.0")
  skip_if_not(soilVAE::vae_available(), "soilVAE Python backend unavailable")
}

test_that("train_pca_vae fits, scores, saves and predicts end to end", {
  skip_on_cran()
  skip_if_no_soilVAE()
  dat <- tiny_spectra(n = 300L, seed = 4L)
  out <- file.path(tempdir(), "e2e_pcavae"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  y <- apply_target_transform(dat$y_raw, "oc")

  expect_no_error(
    train_pca_vae(dat$X, y, family_id = "FAM", prop = "oc", out_dir = out,
                  n_components = 10L, n_seeds = 2L, hidden_enc = 16L,
                  hidden_dec = 16L, latent_dim = 4L, epochs = 3L, patience = 2L,
                  min_n = 20L)
  )
  expect_trained_artefacts(out, "FAM", "oc", "_pcavae.rds")
  base <- file.path(out, "FAM", "models", "oc")
  expect_true(all(file.exists(paste0(base, c("_vae1.keras", "_vae2.keras")))))

  # the saved bundle holds the PCA basis only, never the training scores
  b <- readRDS(paste0(base, "_pcavae.rds"))
  expect_null(b$pca$scores)
  expect_equal(b$n_seeds, 2L)
  expect_true(all(autoSpectra:::.is_publishable(list.files(dirname(base)))))

  # prediction: mean of the members, intervals around it, in-domain check
  r <- autoSpectra:::.predict_one_pca_vae("FAM", "oc", dat$X[1:5, ], out,
                                          with_pi = TRUE, confidence = 0.9)
  expect_length(r$y, 5L)
  expect_true(all(is.finite(r$y)))
  expect_true(all(r$pi_lo < r$y & r$y < r$pi_hi))
  r1 <- autoSpectra:::.predict_one_pca_vae("FAM", "oc", dat$X[1, , drop = FALSE], out)
  expect_equal(r1$y, r$y[1], tolerance = 1e-5)
  info <- jsonlite::read_json(paste0(base, "_metrics.json"))
  expect_equal(info$latent$method, "q_statistic")
  expect_true(is.finite(unlist(info$latent$ucl_q)))
})

test_that("train_sensor_model makes a family that predict_soil() can use", {
  skip_on_cran()
  skip_if_no_soilVAE()
  set.seed(9)
  n <- 260L; wl <- seq(1350, 2550, by = 2)
  shape <- exp(-((wl - 1900) / 120)^2)
  carbon <- stats::runif(n, 0.2, 4)
  refl <- t(vapply(carbon, function(cc) 0.5 - 0.05 * cc * shape + stats::rnorm(length(wl), sd = 0.002),
                   numeric(length(wl))))
  df <- data.frame(Soil_ID = paste0("s", seq_len(n)), refl)
  names(df)[-1] <- wl
  out <- file.path(tempdir(), "e2e_sensor"); unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_error(train_sensor_model(df, carbon, "oc", "OSSL_VisNIR_Cubist", out_dir = out),
               "built-in")
  # few soils: "auto" picks Cubist
  mc <- train_sensor_model(df, carbon, "oc", "SmallSensor", out_dir = out, n_components = 10L)
  expect_identical(get_family("SmallSensor", model_dir = out)$algorithm, "cubist_pca")
  expect_true(is.finite(mc$R2))
  expect_true(all(is.finite(predict_soil(df[1:3, ], "SmallSensor", model_dir = out)$oc)))
  expect_error(train_sensor_model(df, carbon, "oc", "SmallSensor", model = "pca_vae",
                                  out_dir = out), "already uses")

  m <- train_sensor_model(df, carbon, "oc", "TestSensor", model = "pca_vae", out_dir = out,
                          n_components = 10L, n_seeds = 2L)
  expect_true(is.finite(m$R2))
  fam <- get_family("TestSensor", model_dir = out)
  expect_identical(fam$algorithm, "pca_vae")
  expect_equal(range(fam$wavegrid), c(1350, 2500))        # OSSL grid within the sensor range
  p <- predict_soil(df[1:4, ], "TestSensor", model_dir = out, with_pi = TRUE)
  expect_true(all(is.finite(p$oc)) && all(p$oc_lo <= p$oc & p$oc <= p$oc_hi))
  ad <- predict_applicability(df[1:4, ], "TestSensor", "oc", model_dir = out)
  expect_identical(attr(ad, "method"), "q_statistic")
})
