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
