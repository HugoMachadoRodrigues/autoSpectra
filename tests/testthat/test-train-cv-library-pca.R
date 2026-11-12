# The library PCA of the CV: each fold's PCA must never see a test-fold sample,
# and the pooled metrics must be aggregates only.

test_that("library PCA excludes every test-fold sample and pooled metrics are aggregates", {
  skip_if_not_installed("Cubist")
  set.seed(1)
  n_lab <- 300; n_unl <- 200; p <- 40
  X_lab <- matrix(rnorm(n_lab * p), n_lab, p)
  X_unl <- matrix(rnorm(n_unl * p), n_unl, p)
  y     <- X_lab[, 1] * 2 + rnorm(n_lab, sd = 0.3)
  g_lab <- paste0("s", seq_len(n_lab))
  # Unlabelled spectra: half are other scans of labelled samples
  g_unl <- c(sample(g_lab, 100), paste0("u", seq_len(100)))
  X_lib <- rbind(X_lab, X_unl); g_lib <- c(g_lab, g_unl)
  folds <- random_kfold(n_lab, k = 3, seed = 42, groups = g_lab)

  seen <- list()
  orig <- autoSpectra:::.fit_pca
  local_mocked_bindings(.fit_pca = function(X, n_components = 120L, scale = FALSE) {
    seen[[length(seen) + 1L]] <<- X
    orig(X, n_components = 5L, scale = scale)
  })

  fam <- list(algorithm = "cubist_pca")
  res <- cv_train_property(X_lab, y, "FAM", "clay.tot", folds, fam,
                           groups = g_lab, X_lib = X_lib, lib_groups = g_lib,
                           pca_scale = TRUE, committees_grid = 1L, k_inner = 2L)

  expect_length(seen, 3L)
  for (i in 1:3) {
    test_rows <- rbind(X_lab[folds[[i]], , drop = FALSE],
                       X_unl[g_unl %in% g_lab[folds[[i]]], , drop = FALSE])
    key <- function(M) apply(round(M, 10), 1, paste, collapse = ",")
    expect_false(any(key(test_rows) %in% key(seen[[i]])))
    # ...and it is larger than the labelled training rows alone
    expect_gt(nrow(seen[[i]]), n_lab - length(folds[[i]]))
  }

  expect_identical(res$pca_fit, "library minus test fold")
  expect_true(all(c("R2", "RMSE", "cor2", "n") %in% names(res$pooled)))
  expect_equal(res$pooled$n, n_lab)
  # Nothing per sample is returned: every element is short
  lens <- unlist(rapply(res, length, how = "unlist"))
  expect_true(all(lens < n_lab))
  expect_null(res$per_fold[[1]]$.oof)
})

test_that("a scaled PCA reconstructs the spectra (Q-statistic stays valid)", {
  set.seed(2)
  X <- matrix(rnorm(200 * 12), 200, 12) %*% diag(seq(0.1, 5, length.out = 12))
  pca <- autoSpectra:::.fit_pca(X, n_components = 12L, scale = TRUE)
  X_hat <- autoSpectra:::.reconstruct_pca(pca, autoSpectra:::.project_pca(pca, X))
  expect_equal(unname(X_hat), unname(X), tolerance = 1e-8)
})

test_that(".fit_split_pca leaves the held-out samples out of the library", {
  set.seed(3)
  X_lib <- matrix(rnorm(300 * 10), 300, 10); g_lib <- rep(paste0("g", 1:150), 2)
  X_tr  <- X_lib[1:100, ]
  seen <- NULL
  orig <- autoSpectra:::.fit_pca
  local_mocked_bindings(.fit_pca = function(X, n_components = 120L, scale = FALSE) {
    seen <<- X; orig(X, n_components = n_components, scale = scale)
  })
  pca <- autoSpectra:::.fit_split_pca(X_tr, n_components = 5L, scale = TRUE,
                                      X_lib = X_lib, lib_groups = g_lib,
                                      exclude_groups = c("g140", "g150"))
  expect_equal(nrow(seen), 296L)                 # 2 groups x 2 scans left out
  expect_equal(pca$n_library, 296L)
  expect_equal(dim(pca$scores), c(100L, 5L))     # scores = the training rows
  expect_equal(pca$scores, autoSpectra:::.project_pca(pca, X_tr))
})

test_that("train_cubist_pca with a library never fits the PCA on test spectra", {
  skip_if_not_installed("Cubist")
  set.seed(4)
  n <- 300; p <- 20
  X <- matrix(rnorm(n * p), n, p); y <- X[, 1] + rnorm(n, sd = 0.2)
  g <- paste0("s", seq_len(n))
  X_lib <- rbind(X, matrix(rnorm(100 * p), 100, p)); g_lib <- c(g, paste0("u", 1:100))
  seen <- NULL
  orig <- autoSpectra:::.fit_pca
  local_mocked_bindings(.fit_pca = function(X, n_components = 120L, scale = FALSE) {
    seen <<- X; orig(X, n_components = 10L, scale = scale)
  })
  out <- withr::local_tempdir()
  train_cubist_pca(X, y, "FAM", "clay.tot", out_dir = out, groups = g,
                   committees_grid = 1L, k_inner = 2L,
                   X_lib = X_lib, lib_groups = g_lib, pca_scale = TRUE)
  test_rows <- autoSpectra:::.grouped_partition(n, c(0.8, 0.2), g, seed = 42L)[[2]]
  key <- function(M) apply(round(M, 10), 1, paste, collapse = ",")
  expect_false(any(key(X[test_rows, ]) %in% key(seen)))
  expect_equal(nrow(seen), nrow(X_lib) - length(test_rows))
  b <- readRDS(file.path(out, "FAM", "models", "clay.tot_cubist.rds"))
  expect_null(b$pca$scores)
  expect_false(all(b$pca$scale == 1))
})
