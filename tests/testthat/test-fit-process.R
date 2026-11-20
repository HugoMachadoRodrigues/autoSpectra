# The GPU path fits the masked soilVAE in inst/python/fit_masked_vae.py and
# loads its weights into the network built by build_maskedVAE(). The two
# definitions must stay identical.

skip_if_no_keras <- function() {
  skip_if_not_installed("keras3")
  skip_if_not(
    tryCatch({ autoSpectra:::.init_vae_python(); TRUE }, error = function(e) FALSE),
    "keras3 Python backend unavailable"
  )
}

test_that("the Python fit script builds the same network as build_maskedVAE()", {
  skip_on_cran()
  skip_if_no_keras()
  skip_if_not_installed("processx")
  script <- system.file("python", "fit_masked_vae.py", package = "autoSpectra")
  skip_if(!nzchar(script), "fit_masked_vae.py not installed")
  res <- processx::run(reticulate::py_exe(), c(script, "--describe", "40", "8"),
                       env = c("current", KERAS_BACKEND = "tensorflow",
                               TF_CPP_MIN_LOG_LEVEL = "3"),
                       error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  json_line <- grep("^\\[\\[", strsplit(res$stdout, "\n")[[1]], value = TRUE)
  py <- jsonlite::fromJSON(json_line)
  mdl <- build_maskedVAE(40L, latent_dim = 8L)
  r <- t(vapply(mdl$layers, function(l)
    c(class(l)[1], as.character(l$count_params())), character(2)))
  r[, 1] <- sub("^keras\\.src\\.layers\\..*\\.", "", r[, 1])
  r[, 1] <- sub("^.*\\.", "", r[, 1])
  expect_equal(nrow(py), nrow(r))
  expect_equal(py[, 1], r[, 1])
  expect_equal(as.integer(py[, 2]), as.integer(r[, 2]))
})

test_that("a failed GPU fit falls back to the CPU", {
  skip_on_cran()
  skip_if_no_keras()
  skip_if_not_installed("processx")
  set.seed(3)
  wg <- seq(350, 2500, by = 25); n <- 120
  X <- matrix(runif(n * length(wg)), n)
  y <- rowMeans(X[, 1:10]) * 10
  out <- tempfile("fitproc_"); dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  old <- options(autoSpectra.vae_device = "gpu",
                 autoSpectra.python = file.path(out, "no-such-python"))
  on.exit(options(old), add = TRUE)
  expect_message(
    train_maskedVAE(X, y, "OSSL_VisNIR", "clay.tot", wavegrid = wg,
                    sensor_type = "visnir", out_dir = out, epochs = 2L,
                    batch_size = 32L, kl_warmup_epochs = 1L, n_aug_copies = 1L,
                    preprocess = "SNV_MASKED"),
    "fitting on the CPU")
  expect_true(file.exists(file.path(out, "OSSL_VisNIR", "models", "clay.tot.keras")))
})
