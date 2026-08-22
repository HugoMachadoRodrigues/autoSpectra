test_that("resolve_target_transform uses property defaults for aligned OSSL targets", {
  expect_identical(resolve_target_transform("oc", y = c(0, 1, 2), method = "auto"), "log1p")
  expect_identical(resolve_target_transform("cec", y = c(2, 5, 10), method = "auto"), "log1p")
  expect_identical(resolve_target_transform("ph.h2o", y = c(4.5, 6.2, 7.1), method = "auto"), "none")
  expect_identical(resolve_target_transform("clay.tot", y = c(10, 25, 40), method = "auto"), "none")
})

test_that("resolve_target_transform falls back to distribution heuristics", {
  y_skewed <- c(rep(0, 10), 1, 2, 5, 10, 20, 50, 100)
  expect_identical(resolve_target_transform("custom.prop", y = y_skewed, method = "auto"), "log1p")

  y_regular <- c(1.2, 1.5, 1.9, 2.1, 2.4, 2.8, 3.0)
  expect_identical(resolve_target_transform("custom.prop", y = y_regular, method = "auto"), "none")
})

test_that("transform_target rejects unresolved auto transforms", {
  expect_error(transform_target(c(1, 2, 3), method = "auto"), "must be resolved")
  expect_error(inverse_target_transform(c(1, 2, 3), method = "auto"), "must be resolved")
})
