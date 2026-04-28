test_that("q_statistic returns row-wise sum of squared residuals", {
  E <- matrix(c(1, 2, -1, 0, 3, 4), nrow = 2, byrow = TRUE)
  expect_equal(q_statistic(E), c(1 + 4 + 1, 0 + 9 + 16))
})

test_that("q_statistic_ucl is positive and finite for a typical residual matrix", {
  set.seed(10)
  E   <- matrix(rnorm(200, sd = 0.05), nrow = 20)
  ucl <- q_statistic_ucl(E, alpha = 0.01)
  expect_true(is.finite(ucl))
  expect_true(ucl > 0)
})

test_that("mahalanobis_d2 reduces to Euclidean^2 when Sigma is identity", {
  Z  <- matrix(c(1, 2, -1, 0,
                 3, 4,  0, 1), nrow = 2, byrow = TRUE)
  mu <- rep(0, 4)
  Si <- diag(4)
  expect_equal(mahalanobis_d2(Z, mu, Si),
               c(1 + 4 + 1 + 0, 9 + 16 + 0 + 1))
})

test_that("mahalanobis_threshold matches qchisq", {
  expect_equal(mahalanobis_threshold(16, alpha = 0.05),
               stats::qchisq(0.95, df = 16))
})
