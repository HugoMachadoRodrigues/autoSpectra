test_that("drop_empty_spectra removes rows with no finite bands", {
  M <- matrix(
    c(NA, NA,
      1,  2,
      NA, 3),
    nrow = 3,
    byrow = TRUE
  )
  rownames(M) <- c("a", "b", "c")
  colnames(M) <- c("350", "352")

  out <- .drop_empty_spectra(M)

  expect_equal(attr(out, "n_dropped_empty"), 1L)
  expect_equal(rownames(out), c("b", "c"))
  expect_equal(nrow(out), 2L)
})

test_that("ossl_visnir_matrix drops empty spectra rows", {
  visnir_df <- data.frame(
    id.layer_uuid_txt = c("id1", "id2"),
    scan_visnir.350_ref = c(NA, 0.2),
    scan_visnir.352_ref = c(NA, 0.3),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  M <- ossl_visnir_matrix(visnir_df)

  expect_equal(attr(M, "n_dropped_empty"), 1L)
  expect_equal(rownames(M), "id2")
  expect_equal(colnames(M), c("350", "352"))
})

test_that("spectral availability report flags zero-spectrum datasets", {
  visnir_df <- data.frame(
    id.layer_uuid_txt = c("id1", "id2", "id3", "id4"),
    dataset.code_ascii_txt = c("LUCAS.SSL", "KSSL.SSL", "KSSL.SSL", "AFSIS1.SSL"),
    scan_visnir.350_ref = c(0.2, 0.3, NA, NA),
    scan_visnir.352_ref = c(0.4, 0.5, NA, NA),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  M <- ossl_visnir_matrix(visnir_df)
  report <- .ossl_spectral_availability_report(
    spec_raw = visnir_df,
    sensor_type = "visnir",
    valid_ids = rownames(M)
  )

  counts <- report$counts

  expect_equal(report$n_datasets_raw, 3L)
  expect_equal(report$n_datasets_valid, 2L)
  expect_equal(report$n_datasets_empty, 1L)
  expect_equal(counts$n_valid[counts$dataset == "LUCAS.SSL"], 1L)
  expect_equal(counts$n_valid[counts$dataset == "KSSL.SSL"], 1L)
  expect_equal(counts$n_empty[counts$dataset == "KSSL.SSL"], 1L)
  expect_equal(counts$status[counts$dataset == "KSSL.SSL"], "partial")
  expect_equal(counts$status[counts$dataset == "AFSIS1.SSL"], "no_valid_spectra")
  expect_true(counts$possible_upstream_inconsistency[counts$dataset == "AFSIS1.SSL"])
})
