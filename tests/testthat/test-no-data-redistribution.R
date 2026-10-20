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
  pkg <- skip_if_no_source_tree()

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

test_that(".is_publishable allows model artefacts and refuses everything else", {
  ok <- c("oc.keras", "oc_cubist.rds", "oc_cubist_error.rds",
          "oc_scaler.rds", "oc_scaler.json", "oc_metrics.json",
          "metrics_summary.json", "cv10_summary.json")
  expect_true(all(autoSpectra:::.is_publishable(ok)))

  # The shapes a per-sample dump would plausibly take must all be refused,
  # because a denylist publisher is what let the Cubist $data blob through.
  bad <- c("oc_residuals.rds", "train_scores.csv", "ossl_visnir_L0_v1.2.csv.gz",
           "fold_indices.rds", "predictions.parquet", "X_train.qs",
           "spectra.feather", "oc_training_data.json", "notes.txt")
  expect_false(any(autoSpectra:::.is_publishable(bad)))
})

test_that(".is_publishable matches on basename, not on the directory path", {
  # A stray dump inside a directory whose name ends in .keras must not pass
  expect_false(autoSpectra:::.is_publishable("models/oc.keras/train_rows.csv"))
  expect_true(autoSpectra:::.is_publishable("models/OSSL_MIR/models/oc.keras"))
  expect_equal(autoSpectra:::.is_publishable(character(0)), logical(0))
})

test_that(".ossl_gz_ok rejects truncated, empty and non-gzip files", {
  skip_if(!nzchar(Sys.which("gzip")), "gzip binary not available")

  # A valid gzip round-trip is accepted. Use high-entropy content: identical
  # repeated lines would compress below the function's 1 KB sanity floor, which
  # is a legitimate guard for OSSL files (smallest is ~10 MB).
  set.seed(42)
  good <- tempfile(fileext = ".csv.gz"); on.exit(unlink(good), add = TRUE)
  con  <- gzfile(good, "wb")
  writeLines(vapply(seq_len(3000),
                    function(i) paste(runif(4), collapse = ","),
                    character(1)), con)
  close(con)
  expect_gt(file.size(good), 1024)
  expect_true(autoSpectra:::.ossl_gz_ok(good))

  # Truncating it must be detected. R's own gzfile() reads a truncated member
  # up to the cut and reports a clean EOF, so this is the case that matters.
  trunc <- tempfile(fileext = ".csv.gz"); on.exit(unlink(trunc), add = TRUE)
  raw   <- readBin(good, "raw", n = file.size(good))
  writeBin(raw[seq_len(floor(length(raw) * 0.6))], trunc)
  expect_false(autoSpectra:::.ossl_gz_ok(trunc))

  empty <- tempfile(fileext = ".csv.gz"); on.exit(unlink(empty), add = TRUE)
  file.create(empty)
  expect_false(autoSpectra:::.ossl_gz_ok(empty))

  junk <- tempfile(fileext = ".csv.gz"); on.exit(unlink(junk), add = TRUE)
  writeBin(as.raw(rep(0L, 5000)), junk)
  expect_false(autoSpectra:::.ossl_gz_ok(junk))

  expect_false(autoSpectra:::.ossl_gz_ok(tempfile()))
})

test_that("ossl_download targets the official OSSL storage bucket", {
  # Users must pull OSSL from the official source documented at
  # docs.soilspectroscopy.org, never from a mirror we host.
  expect_equal(autoSpectra:::.OSSL_BASE,
               "https://storage.googleapis.com/soilspec4gg-public/")

  # And the cache must be user-local, outside the package tree
  expect_false(grepl("autoSpectra/(R|inst|data)/", ossl_cache_dir()))
})

test_that(".Rbuildignore excludes every directory that can hold data or artefacts", {
  pkg <- skip_if_no_source_tree()
  rbi <- readLines(file.path(pkg, ".Rbuildignore"), warn = FALSE)

  # Directories that hold OSSL-derived output, local scratch, or session state.
  # If any of these ever ships inside the package tarball, the package is
  # distributing something other than code.
  must_exclude <- c("\\.claude", "\\.git", "dist", "models",
                    "backup_models", "results", "data test")
  for (d in must_exclude)
    expect_true(any(grepl(d, rbi, fixed = TRUE)),
                info = paste("missing from .Rbuildignore:", d))
})

test_that("the package source tree carries no data files outside ignored paths", {
  pkg <- skip_if_no_source_tree()
  # Only the directories R CMD build actually copies from
  shipped <- unlist(lapply(c("R", "man", "inst", "data", "src", "tests"),
                           function(d) {
                             p <- file.path(pkg, d)
                             if (!dir.exists(p)) return(character(0))
                             list.files(p, recursive = TRUE)
                           }))
  data_files <- grep("\\.(csv|rda|RData|rds|qs2?|parquet|feather|xlsx|xls|gz|zip)$",
                     shipped, ignore.case = TRUE, value = TRUE)
  expect_equal(data_files, character(0))
})
