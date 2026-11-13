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

test_that("NeoSpectra soils join OSSL lab rows through XS/XN scan ids", {
  scan  <- c("235157XS", "235158XN", "999", NA, "235159XS")
  layer <- c("L1",       "L2",       "L3",  "L4", "L5")
  # OSSL lab rows in another order, plus a layer with no scan
  got <- autoSpectra:::.neospectra_id_from_scan(scan, layer, c("L5", "L3", "L1", "L9", "L2"))
  expect_identical(got, c("235159", NA, "235157", NA, "235158"))
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

test_that("VisNIR grids: 350-2500 nm for the VAE, 400-2500 nm for the OSSL replication", {
  expect_equal(range(get_family("OSSL_VisNIR")$wavegrid),        c(350, 2500))
  expect_equal(range(get_family("OSSL_VisNIR_Cubist")$wavegrid), c(400, 2500))
  expect_equal(range(get_family("OSSL_VisNIR_PCAVAE")$wavegrid), c(400, 2500))
})

test_that("MIR families use the 600-4000 cm-1 grid", {
  expect_equal(range(get_family("OSSL_MIR")$wavegrid),        c(600, 4000))
  expect_equal(range(get_family("OSSL_MIR_Cubist")$wavegrid), c(600, 4000))
})

test_that("Cubist family pipelines match Safanelli (SNV alone)", {
  # OSSL VisNIR models work on reflectance (no absorbance step)
  expect_equal(get_family("OSSL_VisNIR_Cubist")$preprocess, c("SNV"))
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
