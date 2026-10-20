test_that("model_registry exposes the four expected families", {
  expect_setequal(names(model_registry),
                  c("OSSL_VisNIR", "OSSL_MIR",
                    "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist"))
})

test_that("VisNIR families use the Safanelli 350-2500 nm grid", {
  expect_equal(range(get_family("OSSL_VisNIR")$wavegrid),       c(350, 2500))
  expect_equal(range(get_family("OSSL_VisNIR_Cubist")$wavegrid), c(350, 2500))
})

test_that("MIR families use the 600-4000 cm-1 grid", {
  expect_equal(range(get_family("OSSL_MIR")$wavegrid),        c(600, 4000))
  expect_equal(range(get_family("OSSL_MIR_Cubist")$wavegrid), c(600, 4000))
})

test_that("Cubist family pipelines match Safanelli (SNV alone)", {
  expect_equal(get_family("OSSL_VisNIR_Cubist")$preprocess,
               c("ABSORBANCE", "SNV"))
  expect_equal(get_family("OSSL_MIR_Cubist")$preprocess,
               c("SNV"))
})

test_that("VAE family pipelines use SG_SMOOTH + SNV_MASKED", {
  expect_true("SNV_MASKED" %in% get_family("OSSL_VisNIR")$preprocess)
  expect_true("SNV_MASKED" %in% get_family("OSSL_MIR")$preprocess)
})

test_that("ossl_l1_log_transform excludes the five non-skewed properties", {
  excluded <- c("clay.tot", "silt.tot", "sand.tot", "ph.h2o", "ph.cacl2")
  expect_false(any(excluded %in% ossl_l1_log_transform))
  expect_true(all(setdiff(ossl_l1_properties, excluded) %in%
                   ossl_l1_log_transform))
})
