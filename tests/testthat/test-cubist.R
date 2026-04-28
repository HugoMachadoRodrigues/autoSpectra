test_that(".fit_pca returns correct dimensions and uses irlba when available", {
  skip_if_not_installed("irlba")
  set.seed(42)
  # Use n > 5000 to trigger the irlba path
  X <- matrix(rnorm(6000 * 100), 6000, 100)
  pc <- autoSpectra:::.fit_pca(X, n_components = 20)
  expect_equal(nrow(pc$rotation), 100)
  expect_equal(ncol(pc$rotation), 20)
  expect_equal(nrow(pc$scores),   6000)
  expect_equal(ncol(pc$scores),   20)
  expect_equal(length(pc$center), 100)
  expect_true(all(is.finite(pc$var_explained)))
  expect_true(sum(pc$var_explained) > 0 && sum(pc$var_explained) <= 1.01)
})

test_that(".fit_pca fallback (prcomp) works for small n", {
  set.seed(7)
  X <- matrix(rnorm(200 * 50), 200, 50)
  pc <- autoSpectra:::.fit_pca(X, n_components = 10)
  expect_equal(nrow(pc$rotation), 50)
  expect_equal(ncol(pc$rotation), 10)
  # var_explained should sum to <= 1 (it covers the retained variance fraction)
  expect_true(sum(pc$var_explained) <= 1.01)
})

test_that(".project_pca projects new data correctly", {
  set.seed(11)
  X_train <- matrix(rnorm(500 * 80), 500, 80)
  X_test  <- matrix(rnorm(50  * 80), 50,  80)
  pc <- autoSpectra:::.fit_pca(X_train, n_components = 15)
  T_new <- autoSpectra:::.project_pca(pc, X_test)
  expect_equal(dim(T_new), c(50, 15))
  expect_true(all(is.finite(T_new)))
})

test_that(".reconstruct_pca produces finite reconstructions near input", {
  set.seed(3)
  X <- matrix(rnorm(300 * 60), 300, 60)
  pc <- autoSpectra:::.fit_pca(X, n_components = 60)
  X_hat <- autoSpectra:::.reconstruct_pca(pc, pc$scores)
  expect_equal(dim(X_hat), dim(X))
  # With all 60 PCs the reconstruction should be very close
  expect_true(max(abs(X - X_hat)) < 1e-6)
})
