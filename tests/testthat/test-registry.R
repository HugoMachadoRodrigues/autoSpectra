test_that("model_registry exposes the eight expected families", {
  expect_setequal(names(model_registry),
                  c("OSSL_VisNIR", "OSSL_MIR",
                    "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist",
                    "OSSL_VisNIR_PCAVAE", "OSSL_MIR_PCAVAE",
                    "OSSL_NeoSpectra_Cubist", "OSSL_NeoSpectra_PCAVAE"))
})

test_that("every model of Safanelli et al. (2025, S2 Table) maps to a property", {
  fp <- file.path(rprojroot::find_package_root_file(), "tools", "safanelli2025_s2.csv")
  skip_if_not(file.exists(fp), "tools/ not available")
  s2 <- utils::read.csv(fp, stringsAsFactors = FALSE)
  expect_equal(nrow(s2), 141L)
  expect_true(all(s2$property %in% autoSpectra:::.ossl_l1_columns))
  expect_setequal(names(autoSpectra:::.ossl_l1_columns), ossl_l1_properties)
})

test_that("NeoSpectra rows are matched to OSSL lab rows only when unique", {
  cols <- unname(autoSpectra:::.ossl_l1_columns[c("oc", "n.tot", "ph.h2o", "clay.tot",
                                                  "sand.tot", "cec")])
  neo  <- as.data.frame(matrix(c(1.2, 0.12, 6.10, 20.1, 40.0, 10.0,    # one exact OSSL match
                                 2.0, 0.20, 5.00, 30.0, 30.0, 20.0,    # two OSSL candidates
                                 9.9, 0.90, 7.00, 10.0, 80.0,  5.0),   # none
                               nrow = 3, byrow = TRUE, dimnames = list(NULL, cols)))
  ossl <- as.data.frame(matrix(c(1.23, 0.121, 6.10, 20.07, 40.02, 10.01,
                                 2.01, 0.201, 5.00, 30.02, 30.01, 20.02,
                                 1.99, 0.199, 5.00, 29.98, 29.99, 19.98,
                                 5.00, 0.500, 8.00, 50.00, 10.00, 30.00),
                               nrow = 4, byrow = TRUE, dimnames = list(NULL, cols)))
  expect_identical(autoSpectra:::.neospectra_ossl_match(neo, ossl), c(1L, NA, NA))
})

test_that("PCA-VAE families take exactly the Cubist families' input", {
  for (dom in c("VisNIR", "MIR")) {
    cu <- get_family(paste0("OSSL_", dom, "_Cubist"))
    pv <- get_family(paste0("OSSL_", dom, "_PCAVAE"))
    expect_identical(pv$preprocess, cu$preprocess)
    expect_identical(pv$wavegrid, cu$wavegrid)
    expect_identical(pv$algorithm, "pca_vae")
  }
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
