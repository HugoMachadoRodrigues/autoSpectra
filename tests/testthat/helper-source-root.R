# Some tests inspect the SOURCE tree -- .Rbuildignore, the R/ sources, what
# R CMD build would copy. Those files exist in the repository but not in an
# installed package: R CMD check runs tests from <pkg>.Rcheck/tests/testthat,
# where no DESCRIPTION sits above the working directory, and
# rprojroot::find_root() errors. Such tests must skip there, not fail.

source_root <- function() {
  tryCatch(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
           error = function(e) NULL)
}

skip_if_no_source_tree <- function() {
  testthat::skip_if_not_installed("rprojroot")
  root <- source_root()
  if (is.null(root) || !dir.exists(file.path(root, "R")))
    testthat::skip("source tree not available (running against an installed package)")
  root
}
