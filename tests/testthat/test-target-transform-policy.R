# Which OSSL properties are modelled on a log1p scale.
#
# Carried over from the pre-0.4.0 resolve_target_transform() API. That function
# is gone -- 0.4.0 declares the policy as a single vector, ossl_l1_log_transform,
# applied by .select_y() at training time and inverted by expm1() both when
# scoring and when predicting. The POLICY is unchanged and still worth pinning,
# so these tests were rewritten against the current API rather than deleted.
#
# The distinction matters: skewed, non-negative, unbounded concentrations are
# modelled on the log scale; bounded or already-symmetric properties are not.

test_that("skewed concentration properties are modelled on the log1p scale", {
  for (p in c("oc", "n.tot", "c.tot"))
    expect_true(p %in% ossl_l1_log_transform, info = p)
})

test_that("bounded and symmetric properties are left untransformed", {
  # pH is already on a log scale; the texture triplet is a bounded percentage.
  for (p in c("ph.h2o", "ph.cacl2", "clay.tot", "silt.tot", "sand.tot"))
    expect_false(p %in% ossl_l1_log_transform, info = p)
})

test_that("the log-transform policy only names known OSSL L1 properties", {
  expect_true(all(ossl_l1_log_transform %in% ossl_l1_properties))
})

test_that("log1p and expm1 round-trip on the values the policy admits", {
  y <- c(0, 0.5, 1.44, 15.85, 78.45)      # spans the observed oc range
  expect_equal(expm1(log1p(y)), y, tolerance = 1e-12)
})

test_that("metrics are scored on the original scale for log1p properties", {
  # A perfect prediction in log space must score as perfect in original units;
  # this is what the expm1() back-transform in the trainers guarantees, and it
  # is why our numbers are comparable with Safanelli's tables.
  set.seed(3)
  y      <- c(0.1, 1, 5, 20, 60)
  m_perf <- metrics_from_y(expm1(log1p(y)), expm1(log1p(y)))
  expect_equal(m_perf$RMSE, 0, tolerance = 1e-9)
  expect_equal(m_perf$R2,   1, tolerance = 1e-9)

  # and an error in log space must be reported in original units, not log units
  yhat <- expm1(log1p(y) + 0.1)
  m    <- metrics_from_y(y, yhat)
  expect_gt(m$RMSE, 0.1)          # far larger than the 0.1 log-scale offset
})
