test_that("apply_snv centres rows to mean 0 and unit RMS", {
  set.seed(1)
  M  <- matrix(rnorm(50, mean = 5, sd = 2), nrow = 5)
  Ms <- apply_snv(M)
  expect_equal(rowMeans(Ms), rep(0, 5), tolerance = 1e-10)
  expect_equal(sqrt(rowMeans(Ms^2)), rep(1, 5), tolerance = 1e-10)
})

test_that("apply_snv_masked uses only covered bands and zeros uncovered ones", {
  set.seed(2)
  M     <- matrix(rnorm(20, mean = 3), nrow = 2)
  mask  <- matrix(c(rep(1L, 6), rep(0L, 4),
                    rep(0L, 4), rep(1L, 6)), nrow = 2, byrow = TRUE)
  Ms    <- apply_snv_masked(M, mask)
  # Uncovered positions must be zero-filled
  expect_true(all(Ms[1, 7:10] == 0))
  expect_true(all(Ms[2, 1:4] == 0))
  # Covered slice must have mean 0 within numerical tolerance
  expect_equal(mean(Ms[1, 1:6]), 0, tolerance = 1e-8)
  expect_equal(mean(Ms[2, 5:10]), 0, tolerance = 1e-8)
})

test_that("apply_pipeline runs the full Safanelli VisNIR sequence", {
  set.seed(3)
  R <- matrix(runif(40, 0.05, 0.5), nrow = 4)
  out <- apply_pipeline(R, c("ABSORBANCE", "SNV"))
  expect_equal(dim(out), dim(R))
  expect_equal(rowMeans(out), rep(0, 4), tolerance = 1e-10)
})

test_that("mask_from_coverage detects NAs introduced by resample_to_grid", {
  M <- cbind(1:5, NA_real_, 1:5)
  mk <- mask_from_coverage(M)
  expect_equal(mk[, 2], rep(0L, 5))
  expect_true(all(mk[, c(1, 3)] == 1L))
})

test_that("spectral_dropout produces 1 + n_copies stacked rows", {
  set.seed(4)
  X  <- matrix(rnorm(20), nrow = 2, ncol = 10)
  wl <- seq(350, 2500, length.out = 10)
  out <- spectral_dropout(X, wavegrid = wl, sensor_type = "visnir",
                          n_copies = 3, seed = 5)
  expect_equal(nrow(out$X),    2 * (3 + 1))
  expect_equal(nrow(out$mask), nrow(out$X))
  expect_equal(dim(out$X), dim(out$mask))
  # First copy is full coverage by construction
  expect_true(all(out$mask[1:2, ] == 1L))
})

test_that("has_spectral_coverage flags rows with no bands at all", {
  mask <- rbind(
    c(1L, 1L, 1L, 1L),   # full coverage
    c(0L, 0L, 1L, 1L),   # partial -- legitimate, the masked VAE handles it
    c(0L, 0L, 0L, 0L),   # no scan for this instrument -- must be dropped
    c(1L, 0L, 0L, 0L)    # a single band still counts
  )
  expect_equal(has_spectral_coverage(mask), c(TRUE, TRUE, FALSE, TRUE))
})

test_that("an all-NA spectrum is caught BEFORE the zero-fill defeats is.finite", {
  # This is the regression. An OSSL layer with no scan on this instrument is
  # all-NA; the pipeline zero-fills, after which is.finite(rowSums(X)) can no
  # longer tell it apart from real data.
  X <- rbind(
    c(0.3, 0.4, 0.5, 0.6),
    c(NA_real_, NA_real_, NA_real_, NA_real_)
  )
  mask <- mask_from_coverage(X)
  keep <- has_spectral_coverage(mask)
  expect_equal(keep, c(TRUE, FALSE))

  # after zero-filling, the downstream filter is blind to the empty row
  X0 <- X; X0[!is.finite(X0)] <- 0
  expect_equal(is.finite(rowSums(X0)), c(TRUE, TRUE))   # both pass -- the bug
  expect_equal(keep & is.finite(rowSums(X0)), c(TRUE, FALSE))  # guard catches it
})

test_that("empty rows become a constant vector after SNV, hence indistinguishable", {
  X <- rbind(c(0.3, 0.5, 0.4, 0.6), rep(NA_real_, 4), rep(NA_real_, 4))
  X0 <- X; X0[!is.finite(X0)] <- 0
  out <- apply_snv(X0)
  # the two empty rows are identical to each other: identical inputs force
  # identical predictions, which is what caps R2 regardless of the model
  expect_equal(out[2, ], out[3, ])
  expect_equal(length(unique(round(out[2, ], 12))), 1L)
})
