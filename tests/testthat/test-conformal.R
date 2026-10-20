# ── absolute_conformal ────────────────────────────────────────────────────────

test_that("absolute_conformal returns named numeric vector of quantiles", {
  set.seed(1)
  y_true <- rnorm(200)
  y_pred <- y_true + rnorm(200, sd = 0.5)
  out <- absolute_conformal(y_true, y_pred)
  expect_type(out, "double")
  expect_named(out, c("q68", "q90", "q95"))
  expect_true(all(is.finite(out)))
  expect_true(all(out >= 0))
})

test_that("absolute_conformal quantiles are non-decreasing", {
  set.seed(2)
  y_true <- rnorm(300)
  y_pred <- y_true + rnorm(300)
  out <- absolute_conformal(y_true, y_pred, confidence = c(0.5, 0.75, 0.90, 0.95))
  expect_true(all(diff(out) >= 0))
})

test_that("absolute_conformal ignores non-finite residuals", {
  y_true <- c(1, 2, Inf,  3, NaN, 4)
  y_pred <- c(0, 0,   0, NA,   0, 0)
  expect_no_error(absolute_conformal(y_true, y_pred))
})

test_that("absolute_conformal accepts custom confidence levels", {
  y_true <- 1:100
  y_pred <- rep(50, 100)
  out <- absolute_conformal(y_true, y_pred, confidence = c(0.50, 0.95))
  expect_length(out, 2)
  expect_named(out, c("q50", "q95"))
})

# ── train_normalized_conformal + predict_with_pi ──────────────────────────────

test_that("train_normalized_conformal returns expected list structure", {
  skip_if_not_installed("Cubist")
  set.seed(42)
  n <- 250; p <- 60; k <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rowSums(X[, 1:5]) + rnorm(n, sd = 0.5)

  pca <- autoSpectra:::.fit_pca(X, n_components = k)
  T_scores <- pca$scores
  cm <- Cubist::cubist(x = T_scores, y = y, committees = 3)

  result <- train_normalized_conformal(pca, cm, X, y,
                                       confidence = 0.68, committees = 3,
                                       k_cv = 3L, seed = 1L)

  expect_type(result, "list")
  expect_named(result, c("error_model", "alpha_q", "confidence",
                         "cv_folds", "conformal_subset_n"))
  expect_true(is.numeric(result$alpha_q) && length(result$alpha_q) == 1)
  expect_true(is.finite(result$alpha_q) && result$alpha_q > 0)
  expect_equal(result$confidence, 0.68)
  expect_equal(result$cv_folds,   3L)
})

test_that("predict_with_pi returns data.frame with correct columns and ordering", {
  skip_if_not_installed("Cubist")
  set.seed(7)
  n <- 250; p <- 60; k <- 10
  X_tr <- matrix(rnorm(n * p), n, p)
  y_tr <- rowSums(X_tr[, 1:5]) + rnorm(n, sd = 0.5)

  pca  <- autoSpectra:::.fit_pca(X_tr, n_components = k)
  T_tr <- pca$scores
  cm   <- Cubist::cubist(x = T_tr, y = y_tr, committees = 3)
  conf <- train_normalized_conformal(pca, cm, X_tr, y_tr,
                                     confidence = 0.68, committees = 3,
                                     k_cv = 3L, seed = 2L)

  X_new <- matrix(rnorm(30 * p), 30, p)
  out   <- predict_with_pi(pca, cm,
                           conf$error_model, conf$alpha_q,
                           X_new)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("y_pred", "delta_hat", "pi_lo", "pi_hi"))
  expect_equal(nrow(out), 30)
  expect_true(all(is.finite(out$y_pred)))
  expect_true(all(out$delta_hat >= 0))
  # Prediction interval must bracket the point prediction
  expect_true(all(out$pi_lo <= out$y_pred))
  expect_true(all(out$y_pred <= out$pi_hi))
})

test_that("predict_with_pi interval width scales with alpha_q", {
  skip_if_not_installed("Cubist")
  set.seed(99)
  n <- 250; p <- 60; k <- 10
  X_tr <- matrix(rnorm(n * p), n, p)
  y_tr <- rowSums(X_tr[, 1:5]) + rnorm(n, sd = 0.5)

  pca  <- autoSpectra:::.fit_pca(X_tr, n_components = k)
  T_tr <- pca$scores
  cm   <- Cubist::cubist(x = T_tr, y = y_tr, committees = 3)

  # Train at two confidence levels to get different alpha_q values
  conf68 <- train_normalized_conformal(pca, cm, X_tr, y_tr,
                                       confidence = 0.68, committees = 3,
                                       k_cv = 3L, seed = 3L)
  conf95 <- train_normalized_conformal(pca, cm, X_tr, y_tr,
                                       confidence = 0.95, committees = 3,
                                       k_cv = 3L, seed = 3L)

  X_new <- matrix(rnorm(20 * p), 20, p)
  pi68 <- predict_with_pi(pca, cm, conf68$error_model, conf68$alpha_q, X_new)
  pi95 <- predict_with_pi(pca, cm, conf95$error_model, conf95$alpha_q, X_new)

  width68 <- pi68$pi_hi - pi68$pi_lo
  width95 <- pi95$pi_hi - pi95$pi_lo

  # 95% intervals must be at least as wide as 68% intervals on average
  expect_true(mean(width95) >= mean(width68))
})
