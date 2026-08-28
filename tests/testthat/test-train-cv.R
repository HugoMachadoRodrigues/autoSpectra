# Cross-validation fold construction and metric aggregation.
#
# Rewritten for 0.4.0. The pre-0.4.0 file tested resolve_lodo_fold_count(),
# which is gone: fold construction is now random_kfold() for the default outer
# loop and stratified_kfold() for leave-one-dataset-out, and aggregation is
# .aggregate_fold_metrics(). The intent -- folds must partition the data, and
# aggregation must summarise folds honestly -- is unchanged.

test_that("random_kfold partitions every sample exactly once", {
  n     <- 97L
  folds <- random_kfold(n, k = 10L, seed = 42L)

  expect_length(folds, 10L)
  expect_setequal(sort(unlist(folds)), seq_len(n))
  expect_equal(length(unlist(folds)), n)          # no sample in two folds
  expect_true(all(lengths(folds) > 0L))
  # balanced to within one sample
  expect_lte(max(lengths(folds)) - min(lengths(folds)), 1L)
})

test_that("random_kfold is reproducible for a given seed", {
  expect_equal(random_kfold(50L, k = 5L, seed = 7L),
               random_kfold(50L, k = 5L, seed = 7L))
  expect_false(identical(random_kfold(50L, k = 5L, seed = 7L),
                         random_kfold(50L, k = 5L, seed = 8L)))
})

test_that("stratified_kfold keeps each dataset within a single fold (LODO)", {
  ds    <- rep(c("KSSL.SSL", "LUCAS.SSL", "AFSIS1.SSL", "CAF.SSL"),
               times = c(40, 30, 20, 10))
  folds <- stratified_kfold(ds, k = 4L, seed = 42L)

  expect_setequal(sort(unlist(folds)), seq_along(ds))
  # the point of LODO: a fold never splits a dataset across the train/test line
  for (f in folds)
    expect_equal(length(unique(ds[f])), length(unique(ds[f])), info = "fold")
  for (d in unique(ds)) {
    holding <- vapply(folds, function(f) any(ds[f] == d), logical(1))
    expect_equal(sum(holding), 1L, info = d)
  }
})

test_that(".aggregate_fold_metrics reports mean and SD across folds", {
  folds <- list(
    list(RMSE = 1.0, Bias = 0.0, R2 = 0.80, RPIQ = 2.0, CCC = 0.90, n_test = 10L),
    list(RMSE = 2.0, Bias = 0.2, R2 = 0.60, RPIQ = 1.0, CCC = 0.70, n_test = 20L)
  )
  agg <- autoSpectra:::.aggregate_fold_metrics(folds)

  expect_equal(agg$RMSE_mean, 1.5)
  expect_equal(agg$RMSE_sd,   stats::sd(c(1, 2)))
  expect_equal(agg$R2_mean,   0.70)
  expect_equal(agg$n_folds,   2L)
  expect_equal(agg$n_test,    30L)
})

test_that(".aggregate_fold_metrics tolerates a fold with missing metrics", {
  folds <- list(
    list(RMSE = 1.0, Bias = 0, R2 = 0.8, RPIQ = 2, CCC = 0.9, n_test = 10L),
    list(RMSE = NA_real_, Bias = NA_real_, R2 = NA_real_,
         RPIQ = NA_real_, CCC = NA_real_, n_test = 5L)
  )
  agg <- autoSpectra:::.aggregate_fold_metrics(folds)
  expect_equal(agg$RMSE_mean, 1.0)     # NA folds dropped, not propagated
  expect_equal(agg$n_folds, 2L)
})
