test_that("resolve_lodo_fold_count expands quick mode to full LODO for small subsets", {
  ds <- c("D1", "D1", "D2", "D3", "D4")

  expect_equal(
    autoSpectra:::resolve_lodo_fold_count(ds, requested_k = 3L, quick = TRUE),
    4L
  )
})

test_that("resolve_lodo_fold_count keeps requested quick folds for larger dataset counts", {
  ds <- rep(paste0("D", 1:6), each = 2)

  expect_equal(
    autoSpectra:::resolve_lodo_fold_count(ds, requested_k = 3L, quick = TRUE),
    3L
  )
})

test_that("resolve_lodo_fold_count respects classic LODO outside quick mode", {
  ds <- c("D1", "D1", "D2", "D3", "D4")

  expect_equal(
    autoSpectra:::resolve_lodo_fold_count(ds, requested_k = 10L, quick = FALSE),
    4L
  )
})

test_that("cv aggregation keeps per-fold dataset labels and repeat summaries", {
  run1 <- list(
    summary = list(
      RMSE_mean = 1.5, RMSE_sd = 0.5,
      Bias_mean = 0.1, Bias_sd = 0.1,
      R2_mean = 0.4, R2_sd = 0.2,
      RPIQ_mean = 1.2, RPIQ_sd = 0.3,
      CCC_mean = 0.5, CCC_sd = 0.2,
      n_folds = 2L
    ),
    repeat_summary = data.frame(
      repeat_id = 1L, fold_seed = 42L, train_seed = 42L, n_folds = 2L,
      RMSE_mean = 1.5, Bias_mean = 0.1, R2_mean = 0.4, RPIQ_mean = 1.2, CCC_mean = 0.5
    ),
    fold_summary = data.frame(
      repeat_id = c(1L, 1L),
      fold = c(1L, 2L),
      fold_seed = c(42L, 42L),
      train_seed = c(42L, 42L),
      held_out_dataset = c("ICRAF.ISRIC", "KSSL.SSL"),
      held_out_sensor = c("ASD FieldSpec Pro FR", "ASD Labspec 2500"),
      train_sensors = c("ASD Labspec 2500", "ASD FieldSpec Pro FR"),
      sensor_unseen = c(TRUE, FALSE),
      n_test = c(10L, 12L),
      RMSE = c(1.0, 2.0),
      Bias = c(0.0, 0.2),
      R2 = c(0.6, 0.2),
      RPIQ = c(1.5, 0.9),
      CCC = c(0.7, 0.3)
    )
  )
  run2 <- list(
    summary = list(
      RMSE_mean = 1.7, RMSE_sd = 0.4,
      Bias_mean = 0.0, Bias_sd = 0.2,
      R2_mean = 0.3, R2_sd = 0.1,
      RPIQ_mean = 1.1, RPIQ_sd = 0.2,
      CCC_mean = 0.4, CCC_sd = 0.1,
      n_folds = 2L
    ),
    repeat_summary = data.frame(
      repeat_id = 2L, fold_seed = 1042L, train_seed = 1042L, n_folds = 2L,
      RMSE_mean = 1.7, Bias_mean = 0.0, R2_mean = 0.3, RPIQ_mean = 1.1, CCC_mean = 0.4
    ),
    fold_summary = data.frame(
      repeat_id = c(2L, 2L),
      fold = c(1L, 2L),
      fold_seed = c(1042L, 1042L),
      train_seed = c(1042L, 1042L),
      held_out_dataset = c("ICRAF.ISRIC", "KSSL.SSL"),
      held_out_sensor = c("ASD FieldSpec Pro FR", "ASD Labspec 2500"),
      train_sensors = c("ASD Labspec 2500", "ASD FieldSpec Pro FR"),
      sensor_unseen = c(TRUE, FALSE),
      n_test = c(10L, 12L),
      RMSE = c(1.3, 2.1),
      Bias = c(-0.1, 0.1),
      R2 = c(0.5, 0.1),
      RPIQ = c(1.4, 0.8),
      CCC = c(0.6, 0.2)
    )
  )

  agg <- autoSpectra:::.aggregate_cv_runs(list(run1, run2))

  expect_equal(agg$n_repeats, 2L)
  expect_equal(agg$n_folds, 2L)
  expect_equal(nrow(agg$repeat_summary), 2L)
  expect_equal(nrow(agg$fold_summary), 4L)
  expect_equal(agg$fold_dataset, c("ICRAF.ISRIC", "KSSL.SSL", "ICRAF.ISRIC", "KSSL.SSL"))
  expect_equal(agg$fold_n_test, c(10L, 12L, 10L, 12L))
  expect_equal(agg$R2_mean, mean(c(0.4, 0.3)))
  expect_equal(agg$repeat_R2, c(0.4, 0.3))
  expect_equal(agg$fold_summary$sensor_unseen, c(TRUE, FALSE, TRUE, FALSE))

  datasets <- agg$dataset_summary
  expect_equal(nrow(datasets), 2L)
  expect_equal(
    datasets$held_out_dataset,
    c("ICRAF.ISRIC", "KSSL.SSL")
  )
  expect_equal(
    datasets$n_fold_results,
    c(2L, 2L)
  )
  expect_equal(
    datasets$n_eval_rows,
    c(20L, 24L)
  )
  expect_equal(
    datasets$n_holdout_rows,
    c(20L, 24L)
  )
  expect_equal(
    datasets$R2_mean,
    c(mean(c(0.6, 0.5)), mean(c(0.2, 0.1)))
  )
  expect_equal(
    datasets$sensor_unseen,
    c(TRUE, FALSE)
  )

  groups <- agg$fold_group_summary
  expect_true(all(c("held_out_dataset", "sensor_status") %in% groups$grouping))
  expect_equal(
    groups$n_fold_results[groups$group == "with_icraf"],
    2L
  )
  expect_equal(
    groups$n_fold_results[groups$group == "without_icraf"],
    2L
  )
  expect_equal(
    groups$n_fold_results[groups$group == "sensor_unseen"],
    2L
  )
  expect_equal(
    groups$n_fold_results[groups$group == "sensor_seen"],
    2L
  )
})
