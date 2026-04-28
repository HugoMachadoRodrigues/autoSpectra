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
