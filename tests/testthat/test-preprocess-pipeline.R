test_that("preprocess_visnir returns SNV-scaled spectra", {
  M <- rbind(
    seq(0.2, 0.8, length.out = 25),
    seq(0.3, 0.9, length.out = 25)^1.1
  )

  out <- preprocess_visnir(M, m = 5, p = 2)

  expect_equal(dim(out), dim(M))
  expect_true(all(is.finite(out)))
  expect_equal(rowMeans(out), c(0, 0), tolerance = 1e-6)
  expect_equal(sqrt(rowMeans(out^2)), c(1, 1), tolerance = 1e-6)
})

test_that("official OSSL families use SG smooth plus SNV preprocessing", {
  expect_equal(
    get_family("OSSL_VisNIR")$preprocess,
    c("ABSORBANCE", "SG_SMOOTH(11,2)", "SNV")
  )
  expect_equal(
    get_family("OSSL_MIR")$preprocess,
    c("SG_SMOOTH(11,2)", "SNV")
  )
  expect_equal(
    get_family("OSSL_VisNIR_Masked")$preprocess,
    c("ABSORBANCE", "SG_SMOOTH(11,2)")
  )
})

test_that("spectral_dropout defaults to fixed spectral profiles", {
  wavegrid <- seq(400, 2500, by = 50)
  X_abs <- rbind(
    seq(0.2, 1.0, length.out = length(wavegrid)),
    seq(0.3, 1.1, length.out = length(wavegrid))
  )

  aug <- spectral_dropout(X_abs, wavegrid = wavegrid, seed = 42)
  profile_ranges <- paste0(round(aug$profiles$wl_lo), "-", round(aug$profiles$wl_hi))

  expect_equal(aug$mode, "profiles")
  expect_equal(aug$n_copies, 3L)
  expect_equal(nrow(aug$profiles), 2L)
  expect_equal(profile_ranges, c("1000-2500", "1350-2500"))
  expect_equal(nrow(aug$X), 3L * nrow(X_abs))
  expect_true(all(aug$mask[(nrow(X_abs) + 1L):(2L * nrow(X_abs)), wavegrid < 1000] == 0))
  expect_true(all(aug$mask[(2L * nrow(X_abs) + 1L):(3L * nrow(X_abs)), wavegrid < 1350] == 0))
})

test_that("spectral_dropout random mode preserves legacy random masking", {
  wavegrid <- seq(400, 2500, by = 50)
  X_abs <- matrix(seq(0.2, 1.0, length.out = length(wavegrid)), nrow = 1)

  aug <- spectral_dropout(
    X_abs,
    wavegrid = wavegrid,
    mode = "random",
    n_aug = 2L,
    wl_lo_min = 900,
    wl_lo_max = 901,
    seed = 42
  )

  expect_equal(aug$mode, "random")
  expect_equal(aug$n_copies, 3L)
  expect_equal(length(aug$cutoffs), 2L)
  expect_true(all(aug$cutoffs >= 900 & aug$cutoffs <= 901))
  expect_equal(nrow(aug$profiles), 2L)
  expect_equal(nrow(aug$X), 3L)
})
