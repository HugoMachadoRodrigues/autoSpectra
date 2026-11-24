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

test_that("hotelling_t2 divides each squared score by its component variance", {
  S <- matrix(c(1, 2,
                3, 4), nrow = 2, byrow = TRUE)
  expect_equal(hotelling_t2(S, sdev = c(1, 2)), c(1 + 4 / 4, 9 + 16 / 4))
  expect_error(hotelling_t2(S, sdev = 1), "one value per column")
})

test_that("hotelling_t2_ucl matches the F formula and tends to chi-squared", {
  n <- 200; k <- 5
  expect_equal(hotelling_t2_ucl(n, k, alpha = 0.01),
               k * (n - 1) * (n + 1) / (n * (n - k)) * stats::qf(0.99, k, n - k))
  expect_equal(hotelling_t2_ucl(1e7, 10), stats::qchisq(0.99, 10), tolerance = 1e-3)
  expect_true(is.na(hotelling_t2_ucl(5, 10)))
})

test_that(".fit_pca records component SDs, and T2 of its own rows averages k", {
  set.seed(3)
  X   <- matrix(rnorm(400 * 30), 400) %*% matrix(rnorm(30 * 30), 30)
  pca <- autoSpectra:::.fit_pca(X, n_components = 6)
  expect_length(pca$sdev, 6L)
  expect_equal(pca$sdev, apply(pca$scores, 2, stats::sd), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(mean(hotelling_t2(pca$scores, pca$sdev)), 6 * (399 / 400), tolerance = 1e-6)
})

test_that(".pca_sdev recovers the SDs of a scaled PCA without a stored sdev", {
  set.seed(4)
  X   <- matrix(rnorm(300 * 20), 300) %*% matrix(rnorm(20 * 20), 20)
  pca <- autoSpectra:::.fit_pca(X, n_components = 20, scale = TRUE)
  old <- pca; old$sdev <- NULL          # a model saved before sdev was recorded
  expect_equal(autoSpectra:::.pca_sdev(old), pca$sdev, tolerance = 1e-8)
  uns <- autoSpectra:::.fit_pca(X, n_components = 5); uns$sdev <- NULL
  expect_null(autoSpectra:::.pca_sdev(uns))  # unscaled: not recoverable
})

test_that("predict_hotelling_t2 flags spectra far along the components", {
  set.seed(5)
  X   <- matrix(rnorm(500 * 15), 500) %*% matrix(rnorm(15 * 15), 15)
  pca <- autoSpectra:::.fit_pca(X, n_components = 4)
  ucl <- stats::quantile(hotelling_t2(pca$scores, pca$sdev), 0.99, names = FALSE)
  far <- matrix(colMeans(X), 1) + 10 * pca$sdev[1] * t(pca$rotation[, 1])
  d   <- predict_hotelling_t2(pca, rbind(X[1:3, ], far), ucl)
  expect_named(d, c("t2", "ucl_t2", "in_domain"))
  expect_false(d$in_domain[4])
  expect_true(d$t2[4] > 90)
  # inside the PCA subspace, so the Q-statistic alone would not flag it
  expect_true(predict_q_statistic(pca, far, 1e-6)$q_stat < 1e-6 * sum(X^2))
})

test_that(".t2_limit reads the metrics, and says once when the limit is missing", {
  expect_equal(autoSpectra:::.t2_limit("F", "oc", list(latent = list(ucl_t2 = 123.4))), 123.4)
  note <- autoSpectra:::.t2_note
  old  <- note$shown; on.exit(note$shown <- old)
  note$shown <- NULL
  expect_message(v <- autoSpectra:::.t2_limit("F", "oc", list(latent = list())),
                 "Q-statistic only")
  expect_true(is.na(v))
  expect_silent(autoSpectra:::.t2_limit("F", "oc", list(latent = list())))
})
