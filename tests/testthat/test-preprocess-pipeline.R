# Preprocessing contracts for the four official families.
#
# Rewritten for 0.4.0. The pre-0.4.0 version asserted a single SG-smooth + SNV
# pipeline shared by every family and referenced OSSL_VisNIR_Masked, which no
# longer exists. 0.4.0 deliberately gives the two algorithms different
# pipelines: the Cubist families reproduce Safanelli exactly (SNV, no
# derivative), while the masked soilVAE families add SG smoothing and compute
# SNV only over the covered bands so a partial-coverage instrument is scaled by
# its own bands rather than by absent ones.

test_that("preprocess_visnir returns SNV-scaled spectra", {
  M <- rbind(
    seq(0.2, 0.8, length.out = 25),
    seq(0.3, 0.9, length.out = 25)^1.1
  )
  out <- preprocess_visnir(M)

  expect_equal(dim(out), dim(M))
  expect_true(all(is.finite(out)))
  expect_equal(rowMeans(out), c(0, 0), tolerance = 1e-6)
  expect_equal(sqrt(rowMeans(out^2)), c(1, 1), tolerance = 1e-6)
})

test_that("the Cubist families use Safanelli's SNV pipeline, with no derivative", {
  expect_equal(get_family("OSSL_VisNIR_Cubist")$preprocess,
               c("ABSORBANCE", "SNV"))
  expect_equal(get_family("OSSL_MIR_Cubist")$preprocess, "SNV")

  for (fid in c("OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist"))
    expect_false(any(grepl("DERIV", get_family(fid)$preprocess)), info = fid)
})

test_that("the masked soilVAE families smooth and use per-mask SNV", {
  expect_equal(get_family("OSSL_VisNIR")$preprocess,
               c("ABSORBANCE", "SG_SMOOTH(11,2)", "SNV_MASKED"))
  expect_equal(get_family("OSSL_MIR")$preprocess,
               c("SG_SMOOTH(11,2)", "SNV_MASKED"))
})

test_that("only VisNIR families convert reflectance to absorbance", {
  # OSSL ships VisNIR as reflectance and MIR already as absorbance.
  expect_true("ABSORBANCE" %in% get_family("OSSL_VisNIR")$preprocess)
  expect_true("ABSORBANCE" %in% get_family("OSSL_VisNIR_Cubist")$preprocess)
  expect_false("ABSORBANCE" %in% get_family("OSSL_MIR")$preprocess)
  expect_false("ABSORBANCE" %in% get_family("OSSL_MIR_Cubist")$preprocess)
})

test_that("every family's pipeline is executable by apply_pipeline", {
  for (fid in c("OSSL_VisNIR", "OSSL_MIR",
                "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist")) {
    fam  <- get_family(fid)
    d    <- length(fam$wavegrid)
    X    <- matrix(runif(3 * d, 0.1, 0.9), nrow = 3)
    mask <- matrix(1L, 3, d)
    out  <- apply_pipeline(X, fam$preprocess, mask = mask)
    expect_equal(dim(out), c(3L, d), info = fid)
    expect_true(all(is.finite(out)), info = fid)
  }
})
