# Publishing every family: KSSL families, zip parts under GitHub's 100 MB
# limit, and downloads that reassemble them.

test_that("KSSL families mirror their OSSL pipeline on the KSSL rows", {
  for (alg in c("VisNIR_Cubist", "MIR_Cubist", "VisNIR_PCAVAE", "MIR_PCAVAE")) {
    k <- get_family(paste0("KSSL_", alg)); o <- get_family(paste0("OSSL_", alg))
    expect_identical(k$id, paste0("KSSL_", alg))
    expect_identical(k$subset, "kssl")
    expect_null(o$subset)
    expect_identical(k[c("algorithm", "wavegrid", "preprocess", "sensor_type")],
                     o[c("algorithm", "wavegrid", "preprocess", "sensor_type")])
  }
  expect_true(all(autoSpectra:::.PUBLISHED_FAMILIES %in% names(model_registry)))
})

test_that("window families are published, the masked soilVAE is not", {
  pub <- autoSpectra:::.PUBLISHED_FAMILIES
  expect_true(all(unlist(autoSpectra:::.AGNOSTIC_WINDOWS) %in% pub))
  expect_false(any(c("OSSL_VisNIR", "OSSL_MIR") %in% pub))
  expect_true(all(autoSpectra:::.DEFAULT_DOWNLOAD %in% pub))
})

test_that("format_predictions keeps the model_family text column", {
  p <- data.frame(Soil_ID = c("a", "b"), oc = c(1.234, 2.345),
                  model_family = c("OSSL_VisNIR_PCAVAE", "OSSL_VisNIR_PCAVAE_800"))
  f <- format_predictions(p)
  expect_identical(f[["Model (window)"]], p$model_family)
  expect_equal(f[[2]], c(1.23, 2.35))
})

test_that("family zip names: exact family, parts in numeric order", {
  keys <- c("x/OSSL_MIR_PCAVAE_models.part10.zip", "x/OSSL_MIR_PCAVAE_models.part2.zip",
            "x/OSSL_MIR_PCAVAE_models.part1.zip", "x/OSSL_MIR_models.zip",
            "x/OSSL_MIR_Cubist_models.zip", "x/KSSL_MIR_PCAVAE_models.part1.zip")
  f <- autoSpectra:::.family_zip_names
  expect_identical(f(keys, "OSSL_MIR"), "x/OSSL_MIR_models.zip")
  expect_identical(basename(f(keys, "OSSL_MIR_PCAVAE")),
                   paste0("OSSL_MIR_PCAVAE_models.part", c(1, 2, 10), ".zip"))
  expect_length(f(keys, "OSSL_VisNIR"), 0L)
})

test_that("parts respect the size limit and never split a property", {
  d <- file.path(tempfile(), "FAM", "models"); dir.create(d, recursive = TRUE)
  on.exit(unlink(dirname(dirname(d)), recursive = TRUE), add = TRUE)
  props <- c("p.ext", "p.ext.m3", "oc", "k.ext", "k.ext.m3")
  for (p in props) {
    writeBin(raw(400), file.path(d, paste0(p, "_vae1.keras")))
    writeBin(raw(400), file.path(d, paste0(p, "_vae2.keras")))
    writeBin(raw(10), file.path(d, paste0(p, "_metrics.json")))
    writeBin(raw(10), file.path(d, paste0(p, "_scaler.rds")))
  }
  summ <- file.path(dirname(d), "metrics_summary.json"); writeBin(raw(10), summ)
  files <- c(summ, list.files(d, full.names = TRUE))
  parts <- autoSpectra:::.partition_family_files(files, max_bytes = 1000)
  expect_setequal(unlist(parts), files)
  expect_true(all(vapply(parts, function(x) sum(file.size(x)), numeric(1)) <= 1000))
  expect_true(summ %in% parts[[1]])
  # every property's files sit in exactly one part (p.ext vs p.ext.m3)
  for (p in props) {
    own <- paste0(p, c("_vae1.keras", "_vae2.keras", "_metrics.json", "_scaler.rds"))
    where <- vapply(parts, function(x) any(basename(x) %in% own), logical(1))
    expect_equal(sum(where), 1L)
    expect_true(all(own %in% basename(parts[[which(where)]])))
  }
  expect_length(autoSpectra:::.partition_family_files(files, max_bytes = 1e9), 1L)
})

test_that("download_ossl_models() reassembles a family published in parts", {
  src <- tempfile(); top <- file.path(src, "Owner-autoSpectra-models-abc1234")
  dir.create(top, recursive = TRUE)
  on.exit(unlink(src, recursive = TRUE), add = TRUE)
  fid <- "KSSL_MIR_PCAVAE"
  for (k in 1:2) {
    stage <- tempfile(); d <- file.path(stage, fid, "models")
    dir.create(d, recursive = TRUE)
    p <- c("oc", "clay.tot")[k]
    file.create(file.path(d, paste0(p, c("_pcavae.rds", "_metrics.json"))))
    old <- setwd(stage)
    utils::zip(file.path(top, sprintf("%s_models.part%d.zip", fid, k)), fid, flags = "-rq")
    setwd(old)
  }
  archive <- tempfile(fileext = ".zip")
  old <- setwd(src); utils::zip(archive, basename(top), flags = "-rq"); setwd(old)
  local_mocked_bindings(
    .zenodo_record_files = function(record)
      list(id = 1L, version = "0.5.0", keys = "Owner/autoSpectra-models-v0.5.0.zip",
           urls = archive),
    .download_to = function(url, dest, timeout_sec, record) file.copy(url, dest, overwrite = TRUE)
  )
  out <- tempfile()
  suppressMessages(download_ossl_models(fid, model_dir = out))
  expect_true(all(file.exists(file.path(out, fid, "models",
                                        c("oc_pcavae.rds", "clay.tot_pcavae.rds")))))
  expect_true(models_available(fid, model_dir = out))
  expect_error(download_ossl_models("NOT_A_FAMILY", model_dir = out))
})
