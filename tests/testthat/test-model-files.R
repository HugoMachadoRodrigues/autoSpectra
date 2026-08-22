test_that("model path resolution prefers .keras and falls back to .h5", {
  base <- tempfile("autospectra-model-")
  on.exit(unlink(paste0(base, c(".keras", ".h5"))), add = TRUE)

  file.create(paste0(base, ".h5"))
  expect_match(.resolve_model_path(base), "\\.h5$")

  file.create(paste0(base, ".keras"))
  expect_match(.resolve_model_path(base), "\\.keras$")
})

test_that("model file listing counts both .keras and .h5", {
  dir <- tempfile("autospectra-model-dir-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  file.create(file.path(dir, "a.keras"))
  file.create(file.path(dir, "b.h5"))
  file.create(file.path(dir, "c.txt"))

  files <- .list_model_files(dir)
  expect_length(files, 2L)
  expect_true(all(grepl("\\.(keras|h5)$", files)))
})
