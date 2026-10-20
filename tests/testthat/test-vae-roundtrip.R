# Save -> overwrite -> load -> predict for the masked soilVAE.
#
# Both halves of this have bitten us:
#
#  * keras3::save_model() refuses to write over an existing file unless
#    overwrite = TRUE. Without it a re-run trains a property for its full
#    epoch budget and then throws the result away with
#    "File '<prop>.keras' already exists". The failure only surfaces after
#    the expensive part, which is why it needs a cheap test.
#
#  * The model embeds a custom `Sampling` layer and a custom `masked_mse`
#    loss, so a bare keras3::load_model() cannot deserialize it. Only
#    load_soilVAE(), which passes custom_objects, can read these files back.
#    A model we can write but not read is useless to publish.

skip_if_no_keras <- function() {
  skip_if_not_installed("keras3")
  skip_if_not(
    tryCatch({ autoSpectra:::.init_vae_python(); TRUE }, error = function(e) FALSE),
    "keras3 Python backend unavailable"
  )
}

test_that("a masked VAE survives save -> overwrite -> load -> predict", {
  skip_if_no_keras()

  d   <- 40L
  mdl <- build_maskedVAE(d, latent_dim = 4L)

  root <- file.path(tempdir(), "vae_roundtrip")
  dir.create(file.path(root, "FAM", "models"), recursive = TRUE,
             showWarnings = FALSE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  base <- file.path(root, "FAM", "models", "probe")

  keras3::save_model(mdl, paste0(base, ".keras"), overwrite = TRUE)
  expect_true(file.exists(paste0(base, ".keras")))

  # The regression: writing again over an existing file must not error
  expect_no_error(
    keras3::save_model(mdl, paste0(base, ".keras"), overwrite = TRUE)
  )

  save_scaler(1.5, 2.5, path_base = base)

  bundle <- load_soilVAE("FAM", "probe", model_dir = root)
  expect_equal(bundle$scaler$mean, 1.5)
  expect_equal(bundle$scaler$sd,   2.5)

  X    <- matrix(rnorm(5 * d), 5, d)
  X_in <- autoSpectra:::.pack_masked_input(X, matrix(1L, 5, d))
  pred <- autoSpectra:::.extract_prediction(bundle$model, X_in) *
            bundle$scaler$sd + bundle$scaler$mean

  expect_length(pred, 5L)
  expect_true(all(is.finite(pred)))
})

test_that("every save_model call in the package passes overwrite = TRUE", {
  # Cheap guard that does not need the Python backend: a new save site added
  # without overwrite would reintroduce the train-then-discard failure.
  pkg <- skip_if_no_source_tree()
  srcs <- list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)
  lines <- unlist(lapply(srcs, readLines, warn = FALSE))

  saves <- grep("save_model\\(", lines, value = TRUE)
  saves <- grep("^\\s*#", saves, value = TRUE, invert = TRUE)

  expect_gt(length(saves), 0L)
  expect_true(all(grepl("overwrite\\s*=\\s*TRUE", saves)),
              info = paste("save_model() without overwrite = TRUE:",
                           paste(trimws(saves[!grepl("overwrite\\s*=\\s*TRUE", saves)]),
                                 collapse = " | ")))
})
