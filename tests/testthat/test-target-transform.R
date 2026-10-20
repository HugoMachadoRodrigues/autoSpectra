# The bounded log1p back-transform.
#
# Under leave-one-dataset-out validation the unbounded expm1() destroyed exactly
# the log1p properties: VisNIR Cubist scored n-weighted R2 of -14.83 on n.tot and
# -1.59 on oc, while the untransformed clay.tot (+0.42) and ph.h2o (+0.45) held
# up. On the same folds median CCC stayed at 0.76-0.80, so the models were not
# predicting noise -- a few log-scale overshoots were being exponentiated into
# enormous absolute errors.

test_that("apply_target_transform follows the property policy", {
  y <- c(0, 1, 5)
  expect_equal(apply_target_transform(y, "oc"), log1p(y))
  expect_equal(apply_target_transform(y, "ph.h2o"), y)
})

test_that("apply_target_transform rejects values outside log1p's domain", {
  y <- c(-2, -1, 0, 3)
  out <- apply_target_transform(y, "oc")
  expect_true(is.na(out[1]))   # log1p(-2) undefined
  expect_true(is.na(out[2]))   # log1p(-1) = -Inf
  expect_equal(out[3:4], log1p(c(0, 3)))
})

test_that("target_bounds pads the observed training range", {
  b <- target_bounds(c(0, 10), pad = 0.05)
  expect_equal(b$lo, -0.5)
  expect_equal(b$hi, 10.5)

  # a degenerate target must not produce NaN bounds
  b0 <- target_bounds(numeric(0))
  expect_true(is.infinite(b0$lo) && is.infinite(b0$hi))
})

test_that("an out-of-domain log-scale overshoot no longer explodes", {
  # Training saw organic carbon up to 78.45% -> log1p ceiling of about 4.37.
  y_train <- log1p(c(0, 0.5, 1.44, 15.85, 78.45))
  bounds  <- target_bounds(y_train)

  # A domain-shifted sample predicted far beyond anything in calibration
  overshoot <- 12

  unbounded <- expm1(overshoot)
  bounded   <- invert_target_transform(overshoot, "oc", bounds = bounds)

  expect_gt(unbounded, 100000)          # ~162,754% carbon: physically absurd
  expect_lt(bounded, 120)               # held near the calibration ceiling
  expect_gt(bounded, 78)                # but not crushed below what was seen
})

test_that("predictions inside the calibration range are untouched", {
  y_train <- log1p(c(0, 1, 10, 50))
  bounds  <- target_bounds(y_train)
  inside  <- log1p(c(0.5, 2, 20))
  expect_equal(invert_target_transform(inside, "oc", bounds = bounds),
               expm1(inside))
})

test_that("back-transformed concentrations are never negative", {
  bounds <- target_bounds(log1p(c(0, 10)))
  expect_gte(min(invert_target_transform(c(-5, -0.9, 0), "oc", bounds = bounds)), 0)
})

test_that("untransformed properties pass through unchanged", {
  bounds <- target_bounds(c(3, 9))
  y <- c(4, 6, 8)
  expect_equal(invert_target_transform(y, "ph.h2o", bounds = bounds), y)
  # but still bounded, so a wild pH cannot be returned
  expect_lte(max(invert_target_transform(c(40), "ph.h2o", bounds = bounds)), 9.5)
})

test_that("bounding rescues R2 when a fold overshoots, as LODO did", {
  # Reconstruct the failure mode: most predictions good, a few overshooting the
  # calibration ceiling. This is what drove n.tot to R2 = -46.8 on one fold.
  set.seed(42)
  n       <- 500
  y_raw   <- c(runif(n - 5, 0, 60), rep(2, 5))
  y_model <- log1p(y_raw)
  bounds  <- target_bounds(y_model[seq_len(n - 5)])

  yhat        <- y_model + rnorm(n, sd = 0.15)
  yhat[(n - 4):n] <- 11                      # five domain-shifted overshoots

  naive <- metrics_from_y(expm1(y_model), expm1(yhat))
  fixed <- metrics_original_scale(y_model, yhat, "oc", bounds = bounds)

  expect_lt(naive$R2, -1)      # the unbounded inverse is catastrophic
  expect_gt(fixed$R2, 0.5)     # bounding restores a usable score
  expect_lt(fixed$RMSE, naive$RMSE / 100)
})

test_that("metrics never clamp the observed values, only the predictions", {
  # Clamping truth would flatter the model by moving the target.
  y_model <- log1p(c(0, 1, 100))          # 100 is outside the bounds below
  bounds  <- target_bounds(log1p(c(0, 1)))
  m <- metrics_original_scale(y_model, y_model, "oc", bounds = bounds)
  # identical inputs, yet not a perfect score, because truth was left alone
  # while the prediction was bounded
  expect_gt(m$RMSE, 0)
})
