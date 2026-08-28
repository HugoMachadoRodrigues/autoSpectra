# Guards against shipping OSSL data inside published model artefacts.
#
# The OSSL database must be obtained by each user from the official source
# (https://storage.googleapis.com/soilspec4gg-public/, see R/ossl_data.R).
# autoSpectra must never redistribute it — not as bundled files, and not
# smuggled inside a model file. Two things have to be stripped before any
# Cubist bundle is serialised:
#
#   1. pca$scores  — 120 PC scores per training sample; .reconstruct_pca()
#                    inverts them back to the spectra at ~99.98% of variance.
#   2. model$data  — Cubist keeps a verbatim CSV blob of the FULL training
#                    matrix (predictors + response) in this slot.

test_that(".strip_pca_scores removes per-sample scores but keeps the basis", {
  set.seed(42)
  X   <- matrix(rnorm(400 * 30), 400, 30)
  pca <- autoSpectra:::.fit_pca(X, n_components = 10)

  expect_false(is.null(pca$scores))          # present before stripping
  out <- autoSpectra:::.strip_pca_scores(pca)

  expect_null(out$scores)
  expect_equal(out$n_train, 400L)            # provenance count retained
  # Projection basis must survive — inference depends on it
  expect_equal(dim(out$rotation), c(30L, 10L))
  expect_equal(length(out$center), 30L)
})

test_that(".strip_pca_scores leaves projection of new data unchanged", {
  set.seed(7)
  X_tr <- matrix(rnorm(500 * 40), 500, 40)
  X_te <- matrix(rnorm(25  * 40), 25,  40)
  pca  <- autoSpectra:::.fit_pca(X_tr, n_components = 12)

  before <- autoSpectra:::.project_pca(pca, X_te)
  after  <- autoSpectra:::.project_pca(autoSpectra:::.strip_pca_scores(pca), X_te)
  expect_equal(before, after)
})

test_that(".strip_cubist_data empties the embedded training blob", {
  skip_if_not_installed("Cubist")
  set.seed(3)
  n <- 500; p <- 8
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("PC", seq_len(p))
  y <- X[, 1] * 2 + rnorm(n, sd = 0.3)
  mdl <- Cubist::cubist(x = X, y = y, committees = 2)

  # Cubist really does embed one text row per training sample
  expect_true(nchar(mdl$data) > 0)
  expect_equal(length(strsplit(mdl$data, "\n")[[1]]), n)

  stripped <- autoSpectra:::.strip_cubist_data(mdl)
  expect_identical(stripped$data, "")
})

test_that(".strip_cubist_data does not change predictions", {
  skip_if_not_installed("Cubist")
  set.seed(11)
  n <- 400; p <- 6
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("PC", seq_len(p))
  y <- X[, 2] * 1.5 - X[, 3] + rnorm(n, sd = 0.2)
  mdl <- Cubist::cubist(x = X, y = y, committees = 3)

  X_new <- matrix(rnorm(30 * p), 30, p); colnames(X_new) <- colnames(X)
  expect_equal(predict(mdl, X_new),
               predict(autoSpectra:::.strip_cubist_data(mdl), X_new))
})

test_that("a serialised Cubist bundle carries no training rows", {
  skip_if_not_installed("Cubist")
  set.seed(5)
  n <- 600; p <- 10
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("PC", seq_len(p))
  y <- X[, 1] + rnorm(n, sd = 0.3)

  bundle <- list(
    model = autoSpectra:::.strip_cubist_data(
      Cubist::cubist(x = X, y = y, committees = 2)),
    pca   = autoSpectra:::.strip_pca_scores(
      autoSpectra:::.fit_pca(X, n_components = 5))
  )

  f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
  saveRDS(bundle, f)
  back <- readRDS(f)

  expect_identical(back$model$data, "")
  expect_null(back$pca$scores)

  # A bundle holding 600x10 doubles would be far larger than this
  expect_lt(file.size(f), 200000)
})

test_that("no OSSL data files are bundled in the package source", {
  pkg <- rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))

  # Anything shipped in the built package lives under data/ or inst/.
  # None of it may be an OSSL dump — users download OSSL themselves via
  # ossl_download(), straight from the official Google Cloud Storage bucket.
  shipped <- unlist(lapply(c("data", "inst"), function(d) {
    p <- file.path(pkg, d)
    if (!dir.exists(p)) return(character(0))
    list.files(p, recursive = TRUE)
  }))

  ossl_dumps <- grep("ossl|soilspec|_L0_v|_L1_v", shipped,
                     ignore.case = TRUE, value = TRUE)
  expect_equal(ossl_dumps, character(0))

  # Nor any bulk data format, whatever it is named
  bulk <- grep("\\.(csv|csv\\.gz|qs2?|parquet|feather|arrow)$", shipped,
               ignore.case = TRUE, value = TRUE)
  expect_equal(bulk, character(0))
})

test_that("ossl_download targets the official OSSL storage bucket", {
  # Users must pull OSSL from the official source documented at
  # docs.soilspectroscopy.org, never from a mirror we host.
  expect_equal(autoSpectra:::.OSSL_BASE,
               "https://storage.googleapis.com/soilspec4gg-public/")

  # And the cache must be user-local, outside the package tree
  expect_false(grepl("autoSpectra/(R|inst|data)/", ossl_cache_dir()))
})
