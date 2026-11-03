# Regression tests for the defects found in the November 2026 audit. Each one
# was invisible to the existing suite because the published metrics were scored
# on full-coverage spectra, rounded values looked plausible, or the failure was
# swallowed by a tryCatch.

wl_visnir <- seq(350, 2500, by = 2)

synthetic_reflectance <- function(n = 6, wl = wl_visnir, seed = 1) {
  set.seed(seed)
  t(vapply(seq_len(n), function(i)
    0.3 + 0.1 * i / n + 0.05 * sin(wl / 180) + stats::rnorm(length(wl), sd = 0.002),
    numeric(length(wl))))
}

# ---- 1. Masked VAE: training copies match prediction input ---------------

test_that("augmented training copies equal what predict_soil() computes", {
  fam <- get_family("OSSL_VisNIR")
  R   <- synthetic_reflectance()
  aug <- spectral_dropout(R, wl_visnir, "visnir", n_copies = 3, seed = 9,
                          steps = fam$preprocess)
  n <- nrow(R)
  for (r in seq_len(nrow(aug$X))) {
    partial <- R[(r - 1L) %% n + 1L, ]
    partial[aug$mask[r, ] == 0L] <- NA          # a real instrument with this range
    expected <- prepare_model_input(matrix(partial, 1L), fam$preprocess)$X[1, ]
    expect_equal(unname(aug$X[r, ]), unname(expected), tolerance = 1e-12)
  }
  # and partial copies really are partial
  expect_true(any(rowSums(aug$mask) < ncol(aug$mask)))
})

# ---- 2. No absorbance spike from uncovered bands -------------------------

test_that("uncovered VisNIR bands do not become an absorbance spike", {
  fam <- get_family("OSSL_VisNIR_Cubist")
  R   <- synthetic_reflectance(n = 3)
  R[, wl_visnir < 400] <- NA                     # a LUCAS / Foss-like start
  X   <- prepare_model_input(R, fam$preprocess)$X
  expect_false(anyNA(X))
  expect_true(all(X[, wl_visnir < 400] == 0))
  # Covered bands keep a unit-scale SNV signal; the old zero-fill-then-log path
  # put +6 at the uncovered bands and squeezed the real signal to ~0.2 SD.
  expect_lt(max(abs(X)), 4)
  expect_equal(apply(X[, wl_visnir >= 400], 1, sd), rep(1, 3), tolerance = 0.01)
})

test_that("a partial spectrum gives the same Cubist input as its covered part alone", {
  fam  <- get_family("OSSL_VisNIR_Cubist")
  R    <- synthetic_reflectance(n = 2)
  full <- prepare_model_input(R, fam$preprocess)$X
  cut  <- R; cut[, wl_visnir < 400] <- NA
  part <- prepare_model_input(cut, fam$preprocess)$X
  ref  <- apply_snv(reflectance_to_absorbance(R[, wl_visnir >= 400]))
  expect_equal(unname(part[, wl_visnir >= 400]), unname(ref), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(part, full)))
})

# ---- 3. Q-statistic limit survives the JSON round trip -------------------

test_that("small thresholds are written to JSON at full precision", {
  base <- tempfile()
  save_metrics_json(base, metrics = list(R2 = 0.9), conf_q = list(),
                    feat_range = data.frame(band = "1", min = 0, max = 1),
                    latent = list(ucl_q = 2.44e-5, tiny = 1.234567891e-9))
  back <- jsonlite::read_json(paste0(base, "_metrics.json"))
  expect_equal(back$latent$ucl_q, 2.44e-5)
  expect_equal(back$latent$tiny, 1.234567891e-9)
})

# ---- 4. confidence is honoured -------------------------------------------

test_that("the conformal quantile grid serves any requested confidence", {
  set.seed(3)
  alpha <- abs(stats::rnorm(5000))
  grid  <- autoSpectra:::.conformal_quantile(alpha, autoSpectra:::.CONFORMAL_LEVELS)
  conf  <- list(levels = as.list(autoSpectra:::.CONFORMAL_LEVELS),
                alpha_grid = as.list(grid), alpha_q = grid[19], confidence = 0.68)
  a68 <- autoSpectra:::.conformal_alpha(conf, 0.68)
  a95 <- autoSpectra:::.conformal_alpha(conf, 0.95)
  expect_equal(a68, stats::qnorm(0.84), tolerance = 0.05)   # |N(0,1)| quantile
  expect_equal(a95, stats::qnorm(0.975), tolerance = 0.05)
  expect_error(autoSpectra:::.conformal_alpha(conf, 0.999), "between")
})

test_that("old single-level models warn instead of silently ignoring confidence", {
  old <- list(alpha_q = 1.1, confidence = 0.68)
  expect_warning(a <- autoSpectra:::.conformal_alpha(old, 0.95), "68%")
  expect_equal(a, 1.1)
  expect_silent(autoSpectra:::.conformal_alpha(old, 0.68))
})

test_that("split-conformal intervals cover close to nominal on held-out data", {
  skip_if_not_installed("Cubist")
  set.seed(11)
  n <- 3000; p <- 30
  X <- matrix(stats::rnorm(n * p), n, p)
  y <- X[, 1] + 0.5 * X[, 2]^2 + stats::rnorm(n, sd = 0.3 + 0.3 * abs(X[, 3]))
  tr <- 1:2000; te <- 2001:3000
  pca <- autoSpectra:::.fit_pca(X[tr, ], n_components = 10)
  mdl <- Cubist::cubist(pca$scores, y[tr], committees = 3)
  conf <- train_normalized_conformal(pca, mdl, X[tr, ], y[tr], committees = 3,
                                     k_cv = 5L, seed = 2L)
  for (lv in c(0.68, 0.90)) {
    pi <- predict_with_pi(pca, mdl, conf$error_model,
                          autoSpectra:::.conformal_alpha(conf, lv), X[te, ])
    cover <- mean(y[te] >= pi$pi_lo & y[te] <= pi$pi_hi)
    expect_equal(cover, lv, tolerance = 0.05)
  }
})

# ---- 5. Units match the OSSL columns -------------------------------------

test_that("every property maps to one OSSL column whose units match its label", {
  expect_setequal(names(autoSpectra:::.ossl_l1_columns), ossl_l1_properties)
  unit_of_col <- function(col) sub("^.*_", "", col)
  label_unit  <- c(w.pct = "%", cmolc.kg = "cmolc/kg", mg.kg = "mg/kg",
                   g.cm3 = "g/cm", ds.m = "dS/m", index = "pH", w.frac = "w.frac")
  for (prop in ossl_l1_properties) {
    u <- unit_of_col(autoSpectra:::.ossl_l1_columns[[prop]])
    expect_true(grepl(label_unit[[u]], ossl_l1_labels[[prop]], fixed = TRUE) ||
                  (u == "w.pct" && grepl("w.pct", ossl_l1_labels[[prop]], fixed = TRUE)),
                info = paste(prop, u, ossl_l1_labels[[prop]]))
  }
})

# ---- 6. Savitzky-Golay actually smooths ----------------------------------

test_that("SG smoothing changes the spectrum and matches prospectr", {
  skip_if_not_installed("prospectr")
  set.seed(5)
  M <- matrix(cumsum(stats::rnorm(400)), 2, 200)
  out <- apply_sg_matrix(M, m = 11, p = 2, d = 0)
  expect_false(isTRUE(all.equal(out, M)))
  ref <- prospectr::savitzkyGolay(M, m = 0, p = 2, w = 11)
  expect_equal(unname(out[, 6:195]), unname(ref), tolerance = 1e-10)
})

test_that("SG keeps NA bands NA and smooths each covered run on its own", {
  set.seed(6)
  M <- matrix(cumsum(stats::rnorm(200)), 1, 200)
  M2 <- M; M2[1, 1:40] <- NA
  out <- apply_sg_matrix(M2, 11, 2, 0)
  expect_identical(is.na(out), is.na(M2))
  alone <- apply_sg_matrix(M[, 41:200, drop = FALSE], 11, 2, 0)
  expect_equal(out[1, 41:200], alone[1, ])
})

# ---- 7. Grouped splits ---------------------------------------------------

test_that("LUCAS revisits of one point share a group", {
  df <- data.frame(Soil_ID = c("a", "b", "c", "d"),
                   dataset.code_ascii_txt = c("LUCAS.SSL", "LUCAS.SSL",
                                              "LUCAS.SSL", "KSSL.SSL"),
                   id.scan_local_c = c("2009.47142646", "2015.47142646",
                                       "2009.11111111", "150132XN"))
  g <- autoSpectra:::.ossl_sample_groups(df)
  expect_equal(g[1], g[2])
  expect_false(g[1] == g[3])
  expect_equal(g[4], "d")
})

test_that("grouped splits never put one group on both sides", {
  g  <- c(rep(1:300, each = 2), 301:700)
  sp <- split_idx(length(g), groups = g)
  expect_length(intersect(g[sp$train], g[sp$test]), 0)
  expect_length(intersect(g[sp$train], g[sp$calib]), 0)
  expect_identical(sort(unname(unlist(sp))), seq_along(g))
  folds <- random_kfold(length(g), k = 5, groups = g)
  owner <- integer(length(g))
  for (k in seq_along(folds)) owner[folds[[k]]] <- k
  expect_true(all(tapply(owner, g, function(o) length(unique(o))) == 1L))
})

# ---- Smaller fixes -------------------------------------------------------

test_that("predictions are clamped to the physical range", {
  expect_equal(invert_target_transform(c(-3, 50, 140), "clay.tot"), c(0, 50, 100))
  expect_equal(invert_target_transform(c(-1, 7, 20), "ph.h2o"), c(0, 7, 14))
  expect_true(is.na(apply_target_transform(256, "silt.tot")))
})

test_that("resampling does not interpolate across a wide gap", {
  wl <- c(seq(400, 1000, by = 2), seq(1400, 2500, by = 2))
  M  <- matrix(1, 1, length(wl))
  out <- resample_to_grid(M, wl, seq(350, 2500, by = 2))
  grid <- seq(350, 2500, by = 2)
  expect_true(all(is.na(out[1, grid > 1000 & grid < 1400])))
  expect_true(all(is.finite(out[1, grid >= 400 & grid <= 1000])))
})

test_that("the model cache key includes the model directory", {
  expect_false(autoSpectra:::.cache_key("a", "F", "oc") ==
                 autoSpectra:::.cache_key("b", "F", "oc"))
})

test_that("masked-VAE inference is deterministic", {
  skip_if_not_installed("keras3")
  skip_if_not(tryCatch({ autoSpectra:::.init_vae_python(); TRUE },
                       error = function(e) FALSE), "keras3 backend unavailable")
  mdl <- build_maskedVAE(20, latent_dim = 4)
  x   <- matrix(stats::rnorm(5 * 40), 5, 40)
  expect_identical(autoSpectra:::.extract_prediction(mdl, x),
                   autoSpectra:::.extract_prediction(mdl, x))
})
