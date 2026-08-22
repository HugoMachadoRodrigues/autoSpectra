test_that("replicate-aware sampling keeps all spectra for selected samples", {
  toy <- data.frame(
    Soil_ID = c("A", "A", "B", "C", "D"),
    `400` = c(0.10, 0.11, 0.20, 0.30, 0.40),
    oc = c(1, 1, 2, 3, 4),
    dataset.code_ascii_txt = c("D1", "D1", "D2", "D3", "D2"),
    scan.visnir.model.name_utf8_txt = c(
      "ASD Labspec 2500",
      "ASD Labspec 2500",
      "Metrohm NIRS XDS RapidContent Analyzer",
      "ASD FieldSpec Pro FR",
      "Metrohm NIRS XDS RapidContent Analyzer"
    ),
    stringsAsFactors = FALSE
  )

  subset_info <- autoSpectra:::select_ossl_training_subset(
    toy,
    sample_n    = 3L,
    sensor_col  = "scan.visnir.model.name_utf8_txt",
    dataset_col = "dataset.code_ascii_txt",
    seed        = 42L
  )

  expect_equal(length(unique(subset_info$data$Soil_ID)), 3L)
  expect_equal(sum(subset_info$data$Soil_ID == "A"), 2L)
  expect_true(isTRUE(subset_info$kept_replicates))
  expect_setequal(
    unique(subset_info$sample_meta$primary_sensor),
    c(
      "ASD Labspec 2500",
      "ASD FieldSpec Pro FR",
      "Metrohm NIRS XDS RapidContent Analyzer"
    )
  )
})

test_that("dataset-balanced sampling spreads the subset across datasets", {
  toy <- data.frame(
    Soil_ID = paste0("S", seq_len(14)),
    `400` = seq_len(14) / 10,
    dataset.code_ascii_txt = c(rep("D1", 10), rep("D2", 2), rep("D3", 2)),
    stringsAsFactors = FALSE
  )

  subset_info <- autoSpectra:::select_ossl_training_subset(
    toy,
    sample_n    = 6L,
    dataset_col = "dataset.code_ascii_txt",
    seed        = 42L
  )

  expect_equal(sum(subset_info$dataset_counts), 6L)
  expect_equal(
    as.integer(subset_info$dataset_counts[c("D1", "D2", "D3")]),
    c(2L, 2L, 2L)
  )
  expect_true(isTRUE(subset_info$all_datasets_kept))
})

test_that("dataset + target-balanced sampling spreads the subset across quantile bins", {
  toy <- data.frame(
    Soil_ID = paste0("S", seq_len(16)),
    `400` = seq_len(16) / 10,
    oc = c(1:8, 11:18),
    dataset.code_ascii_txt = rep(c("D1", "D2"), each = 8),
    stringsAsFactors = FALSE
  )

  subset_info <- autoSpectra:::select_ossl_training_subset(
    toy,
    sample_n    = 8L,
    dataset_col = "dataset.code_ascii_txt",
    target_col  = "oc",
    target_bins = 4L,
    seed        = 42L
  )

  bin_tab <- with(
    subset_info$sample_meta,
    table(primary_dataset, target_bin)
  )

  expect_true(isTRUE(subset_info$used_target_balance))
  expect_equal(sum(subset_info$dataset_counts), 8L)
  expect_equal(as.integer(subset_info$dataset_counts[c("D1", "D2")]), c(4L, 4L))
  expect_equal(bin_tab["D1", c("Q1", "Q2", "Q3", "Q4")], rep(1L, 4L))
  expect_equal(bin_tab["D2", c("Q1", "Q2", "Q3", "Q4")], rep(1L, 4L))
})

test_that("1000-sample VisNIR subset covers all sensors in the cached OSSL data", {
  cache_dir <- autoSpectra::ossl_cache_dir()
  visnir_path <- file.path(cache_dir, "ossl_visnir_L0_v1.2.csv.gz")
  soillab_path <- file.path(cache_dir, "ossl_soillab_L1_v1.2.csv.gz")

  skip_if_not(file.exists(visnir_path), "Cached OSSL VisNIR file not available")
  skip_if_not(file.exists(soillab_path), "Cached OSSL soillab file not available")

  joined <- autoSpectra::ossl_prepare(
    sensor_type = "visnir",
    cache_dir = cache_dir,
    properties = "oc",
    download_if_missing = FALSE
  )

  subset_info <- autoSpectra:::select_ossl_training_subset(
    joined,
    sample_n    = 1000L,
    sensor_col  = "scan.visnir.model.name_utf8_txt",
    dataset_col = "dataset.code_ascii_txt",
    target_col  = "oc",
    seed        = 42L
  )

  normalize_sensor <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | !nzchar(x)] <- "Unknown / unspecified sensor"
    x
  }

  full_sensors <- unique(normalize_sensor(joined$scan.visnir.model.name_utf8_txt))
  sub_sensors  <- unique(normalize_sensor(subset_info$data$scan.visnir.model.name_utf8_txt))

  expect_equal(length(unique(subset_info$data$Soil_ID)), 1000L)
  expect_true(isTRUE(subset_info$all_sensors_kept))
  expect_setequal(sub_sensors, full_sensors)
  expect_true(max(subset_info$dataset_counts) - min(subset_info$dataset_counts) <= 1L)
  expect_true(isTRUE(subset_info$used_target_balance))
  expect_true(all(vapply(
    split(subset_info$sample_meta$target_bin, subset_info$sample_meta$primary_dataset),
    function(x) length(unique(x[x != "missing"])) >= 4L,
    logical(1)
  )))
})
