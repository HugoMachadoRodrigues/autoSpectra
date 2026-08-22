test_that("prepare_winsorized_target clips tails and reports metadata", {
  y <- c(0, 1, 2, 3, 100)
  res <- prepare_winsorized_target(y, lower_prob = 0.2, upper_prob = 0.8)

  expect_true(isTRUE(res$applied))
  expect_equal(length(res$y), length(y))
  expect_true(is.finite(res$lower))
  expect_true(is.finite(res$upper))
  expect_lte(min(res$y, na.rm = TRUE), res$lower)
  expect_gte(max(res$y, na.rm = TRUE), res$upper)
  expect_equal(res$n_lower, 1L)
  expect_equal(res$n_upper, 1L)
})

test_that("prepare_winsorized_target can be disabled cleanly", {
  y <- c(1, 2, 3)
  res <- prepare_winsorized_target(y, lower_prob = NULL, upper_prob = NULL)

  expect_false(isTRUE(res$applied))
  expect_equal(res$y, y)
  expect_equal(res$n_lower, 0L)
  expect_equal(res$n_upper, 0L)
})

test_that("flag_global_spectral_outliers flags the requested upper tail fraction", {
  M <- rbind(
    matrix(0, nrow = 9, ncol = 4),
    c(100, 100, 100, 100)
  )
  out <- flag_global_spectral_outliers(M, frac = 0.1)

  expect_equal(length(out$flag), nrow(M))
  expect_equal(sum(out$flag), 1L)
  expect_true(out$flag[10])
})

test_that("prepare_texture_triplet_filter removes impossible textural rows", {
  df <- data.frame(
    clay.tot = c(30, 0, 40, 20, 10),
    silt.tot = c(30, 0, 30, 40, NA),
    sand.tot = c(40, 0, 80, 40, 90)
  )

  out <- prepare_texture_triplet_filter(df, sum_tol = 5)

  expect_true(isTRUE(out$applied))
  expect_equal(out$n_rows, 5L)
  expect_equal(out$n_valid, 1L)
  expect_equal(out$n_all_zero, 1L)
  expect_equal(out$n_bad_sum, 2L)
  expect_equal(out$n_missing_component, 1L)
  expect_equal(out$mask, c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("is_textural_property recognises only the texture triplet", {
  expect_true(is_textural_property("clay.tot"))
  expect_true(is_textural_property("silt.tot"))
  expect_true(is_textural_property("sand.tot"))
  expect_false(is_textural_property("oc"))
})
