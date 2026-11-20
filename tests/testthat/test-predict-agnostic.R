# predict_soil_agnostic() sends each spectrum to the widest window it covers.

.spec_df <- function(ranges, step = 2) {
  wl <- seq(350, 2500, by = step)
  M <- t(vapply(ranges, function(r) ifelse(wl >= r[1] & wl <= r[2],
                                           runif(length(wl)), NA_real_),
                numeric(length(wl))))
  df <- as.data.frame(M); names(df) <- as.character(wl)
  cbind(Soil_ID = paste0("s", seq_along(ranges)), df, stringsAsFactors = FALSE)
}

test_that("VisNIR spectra are routed to the widest covered window", {
  df <- .spec_df(list(c(350, 2500), c(400, 2500), c(800, 2500),
                      c(950, 2500), c(1350, 2500), c(1400, 2500), c(400, 1800)))
  expect_equal(autoSpectra:::.agnostic_route(df, "visnir"),
               c("OSSL_VisNIR_PCAVAE", "OSSL_VisNIR_PCAVAE",
                 "OSSL_VisNIR_PCAVAE_800", "OSSL_VisNIR_PCAVAE_1000",
                 "OSSL_VisNIR_PCAVAE_1350", NA, NA))
})

test_that("MIR spectra are routed to the full or the 700-3800 window", {
  wl <- seq(600, 4000, by = 2)
  mk <- function(lo, hi) ifelse(wl >= lo & wl <= hi, runif(length(wl)), NA_real_)
  df <- as.data.frame(rbind(mk(600, 4000), mk(650, 4000), mk(600, 3800), mk(800, 4000)))
  names(df) <- as.character(wl); df <- cbind(Soil_ID = paste0("m", 1:4), df)
  expect_equal(autoSpectra:::.agnostic_route(df, "mir"),
               c("OSSL_MIR_PCAVAE", "OSSL_MIR_PCAVAE_700_3800",
                 "OSSL_MIR_PCAVAE_700_3800", NA))
})

test_that("predictions are assembled back in the input order", {
  df <- .spec_df(list(c(1350, 2500), c(400, 2500), c(800, 2500), c(2000, 2500)))
  testthat::local_mocked_bindings(
    predict_soil = function(df, family_id, ...) {
      data.frame(Soil_ID = df$Soil_ID, oc = nchar(family_id) + seq_len(nrow(df)) / 10)
    })
  expect_warning(out <- predict_soil_agnostic(df, "visnir", properties = "oc"),
                 "cover none")
  expect_equal(out$Soil_ID, df$Soil_ID)
  expect_equal(out$model_family, c("OSSL_VisNIR_PCAVAE_1350", "OSSL_VisNIR_PCAVAE",
                                   "OSSL_VisNIR_PCAVAE_800", NA))
  expect_equal(out$oc, c(nchar("OSSL_VisNIR_PCAVAE_1350") + 0.1,
                         nchar("OSSL_VisNIR_PCAVAE") + 0.1,
                         nchar("OSSL_VisNIR_PCAVAE_800") + 0.1, NA))
})
