# audit_release_archive() is the independent check on the bytes we publish.
# The write-side guards and the publishing allowlist are promises made by code
# that could regress; this opens a built archive and verifies the property.
#
# An auditor that has not been red-teamed is only a promise too, so these tests
# plant each known leak and require it to be caught, AND require a clean bundle
# to pass. The false-positive half matters as much: an auditor that fails every
# clean release trains people to ignore it.

make_bundle <- function(strip_data, strip_scores, n = 3000L, p = 30L) {
  set.seed(11)
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("PC", seq_len(p))
  y <- X[, 1] + rnorm(n, 0, 0.3)
  m   <- Cubist::cubist(x = X, y = y, committees = 2)
  pca <- autoSpectra:::.fit_pca(X, n_components = 10)
  list(model = if (strip_data)   autoSpectra:::.strip_cubist_data(m) else m,
       pca   = if (strip_scores) autoSpectra:::.strip_pca_scores(pca) else pca)
}

zip_bundle <- function(bundle) {
  root <- file.path(tempdir(), paste0("ra_", sample.int(1e6, 1)))
  d    <- file.path(root, "FAM", "models")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, file.path(d, "probe_cubist.rds"))
  z <- file.path(root, "archive.zip")
  zip::zip(z, files = "FAM/models/probe_cubist.rds", root = root)
  z
}

test_that("a clean bundle passes the audit", {
  skip_if_not_installed("Cubist"); skip_if_not_installed("zip")
  r <- audit_release_archive(
    zip_bundle(make_bundle(strip_data = TRUE, strip_scores = TRUE)),
    quiet = TRUE)
  expect_true(r$ok)
  expect_equal(r$findings, character(0))
})

test_that("Cubist's embedded training matrix is caught", {
  skip_if_not_installed("Cubist"); skip_if_not_installed("zip")
  r <- audit_release_archive(
    zip_bundle(make_bundle(strip_data = FALSE, strip_scores = TRUE)),
    quiet = TRUE)
  expect_false(r$ok)
  expect_true(any(grepl("\\$data", r$findings)))
})

test_that("per-sample PCA scores are caught", {
  skip_if_not_installed("Cubist"); skip_if_not_installed("zip")
  r <- audit_release_archive(
    zip_bundle(make_bundle(strip_data = TRUE, strip_scores = FALSE)),
    quiet = TRUE)
  expect_false(r$ok)
  expect_true(any(grepl("scores", r$findings)))
})

test_that("a file outside the publishing allowlist is caught", {
  skip_if_not_installed("zip")
  root <- file.path(tempdir(), paste0("ra_deny_", sample.int(1e6, 1)))
  d    <- file.path(root, "FAM", "models")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("Soil_ID,score", "abc,1.0"), file.path(d, "train_rows.csv"))
  z <- file.path(root, "archive.zip")
  zip::zip(z, files = "FAM/models/train_rows.csv", root = root)

  r <- audit_release_archive(z, quiet = TRUE)
  expect_false(r$ok)
  expect_true(any(grepl("allowlist", r$findings)))
})

test_that("verbose rule text is NOT mistaken for tabular data", {
  # The real MIR clay.tot $output is 20,810 lines -- 966 rules over 15
  # committees, each printing a 120-coefficient equation across wrapped
  # continuation lines. Judging character blobs by size flagged every clean
  # bundle, so the auditor looks for delimited numeric records instead.
  rule_text <- paste(
    c(rep("    then", 200),
      rep("\toutcome = 0.93 - 0.133 PC3 + 0.087 PC2 + 0.0134 PC1 - 0.131 PC8", 8000),
      rep("  Rule 7/47: [331 cases, mean 44.5, range 8.3 to 93.2]", 400)),
    collapse = "\n")
  expect_equal(autoSpectra:::.audit_r_object(rule_text, "$output"),
               character(0))

  # while an actual delimited record dump is caught
  rows <- paste(replicate(50, paste(round(rnorm(40), 6), collapse = ",")),
                collapse = "\n")
  expect_true(length(autoSpectra:::.audit_r_object(rows, "$data")) > 0L)
})

test_that("aggregate statistics of band or latent size are not flagged", {
  # feat_range spans the band grid (1076 VisNIR / 1701 MIR) and the latent
  # covariance is 32x32; both are legitimate model metadata.
  ok <- list(feat_range = list(min = runif(1701), max = runif(1701)),
             latent = list(mu = runif(32), Sigma = matrix(runif(32 * 32), 32)))
  expect_equal(autoSpectra:::.audit_r_object(ok, ""), character(0))

  # but a per-sample vector is
  expect_true(
    length(autoSpectra:::.audit_r_object(list(resid = runif(80000)), "")) > 0L)
})
