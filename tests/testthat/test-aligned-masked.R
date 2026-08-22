test_that("target transforms round-trip on the raw scale", {
  y <- c(0, 0.5, 2, 10)

  expect_equal(
    inverse_target_transform(transform_target(y, "none"), "none"),
    y
  )
  expect_equal(
    inverse_target_transform(transform_target(y, "log1p"), "log1p"),
    y,
    tolerance = 1e-8
  )
  expect_equal(
    inverse_target_transform(transform_target(y, "sqrt"), "sqrt"),
    y,
    tolerance = 1e-8
  )
})

test_that("feature_minmax falls back to positional band ids when colnames are missing", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  rng <- feature_minmax(X)

  expect_equal(rng$band, c("1", "2"))
  expect_equal(rng$min, c(1, 3))
  expect_equal(rng$max, c(2, 4))
})

test_that("aligned masked payload keeps one reconstruction target per domain", {
  X <- matrix(seq_len(24), nrow = 6)
  M <- matrix(1, nrow = 6, ncol = 4)
  ds <- rep(c("A", "B"), each = 3)
  y <- c(1, 2, 3, 4, 5, 6)

  payload <- autoSpectra:::.build_aligned_fit_payload(
    X_spec = X,
    mask = M,
    y_z = scale(log1p(y))[, 1],
    dataset_labels = ds,
    y_raw = y,
    n_bins = 3L
  )

  expect_length(payload$x, 3L)
  expect_length(payload$y, 3L)
  expect_equal(dim(payload$y[[1L]]), dim(X))
  expect_equal(dim(payload$y[[2L]]), dim(X))
  expect_length(payload$sample_weight, 3L)
  expect_equal(length(payload$sample_weight[[1L]]), nrow(X))
  expect_equal(length(payload$sample_weight[[3L]]), nrow(X))
})

test_that("dataset comparison reports per-dataset deltas for aligned experiments", {
  ref <- list(
    dataset_summary = data.frame(
      held_out_dataset = c("ICRAF.ISRIC", "KSSL.SSL"),
      R2_mean = c(0.1, 0.2),
      RMSE_mean = c(2, 3),
      Bias_mean = c(0.3, -0.4),
      CCC_mean = c(0.5, 0.6),
      RPIQ_mean = c(0.7, 0.8),
      stringsAsFactors = FALSE
    )
  )
  cand <- list(
    dataset_summary = data.frame(
      held_out_dataset = c("ICRAF.ISRIC", "KSSL.SSL"),
      R2_mean = c(0.4, 0.1),
      RMSE_mean = c(1.5, 3.5),
      Bias_mean = c(0.1, -0.2),
      CCC_mean = c(0.7, 0.5),
      RPIQ_mean = c(0.9, 0.7),
      stringsAsFactors = FALSE
    )
  )

  cmp <- autoSpectra:::compare_cv_dataset_summaries(
    ref,
    cand,
    reference_label = "baseline",
    candidate_label = "candidate"
  )

  expect_equal(nrow(cmp), 2L)
  expect_equal(cmp$held_out_dataset, c("ICRAF.ISRIC", "KSSL.SSL"))
  expect_equal(cmp$R2_mean_delta, c(0.3, -0.1))
  expect_equal(cmp$RMSE_mean_delta, c(-0.5, 0.5))
})

test_that("registry exposes the aligned masked family", {
  fam <- get_family("OSSL_VisNIR_Masked_Aligned")

  expect_true(isTRUE(fam$masked))
  expect_true(isTRUE(fam$aligned_masked))
  expect_equal(fam$compare_family, "OSSL_VisNIR_Masked")
  expect_equal(fam$alignment_method, "conditional_coral")
  expect_equal(fam$preprocess, c("ABSORBANCE", "SG_SMOOTH(11,2)"))
})
