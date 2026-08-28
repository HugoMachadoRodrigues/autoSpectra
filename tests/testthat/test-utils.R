test_that("metrics_from_y returns sensible RMSE / R2 for perfect prediction", {
  y <- 1:10
  m <- metrics_from_y(y, y)
  expect_equal(m$RMSE, 0)
  expect_equal(m$R2, 1)
  expect_equal(m$Bias, 0)
})

test_that("get_wavelengths picks numeric column names only", {
  df <- data.frame(Soil_ID = "A",
                   `350` = 0.1, `400` = 0.2,
                   oc = 1.0, ph.h2o = 5.0,
                   check.names = FALSE)
  wl <- get_wavelengths(df)
  expect_equal(wl$wl, c(350, 400))
  expect_setequal(wl$cols, c("350", "400"))
})

test_that("resample_to_grid linearly interpolates and pads with NA", {
  M  <- matrix(c(0, 10, 20), nrow = 1)
  out <- resample_to_grid(M, src_wl = c(100, 200, 300),
                          target_wl = c(150, 250, 400))
  expect_equal(unname(out[1, 1]), 5)
  expect_equal(unname(out[1, 2]), 15)
  expect_true(is.na(out[1, 3]))
})

test_that("split_idx returns three disjoint sets covering 1:n", {
  sp <- split_idx(100, seed = 1)
  all_idx <- sort(c(sp$train, sp$calib, sp$test))
  expect_equal(all_idx, 1:100)
  expect_length(intersect(sp$train, sp$calib), 0)
  expect_length(intersect(sp$train, sp$test),  0)
  expect_length(intersect(sp$calib, sp$test),  0)
  # 60/20/20 split for n=100
  expect_length(sp$train, 60)
  expect_length(sp$calib, 20)
  expect_length(sp$test,  20)
})

test_that("feature_minmax falls back to positional band ids without colnames", {
  M <- matrix(c(1, 2, 3, 4), nrow = 2)   # no dimnames
  fr <- feature_minmax(M)
  expect_equal(nrow(fr), 2L)
  expect_equal(fr$band, c("1", "2"))
  expect_equal(fr$min, c(1, 3))
  expect_equal(fr$max, c(2, 4))
})

test_that("feature_minmax keeps column names when present", {
  M <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(NULL, c("350", "352")))
  fr <- feature_minmax(M)
  expect_equal(fr$band, c("350", "352"))
})
