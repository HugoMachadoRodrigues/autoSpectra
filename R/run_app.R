# R/run_app.R — Shiny app launcher for autoSpectra

#' Launch the autoSpectra Shiny application
#'
#' Opens the interactive soil spectral prediction interface in the default
#' browser. Models must have been trained and saved to \code{model_dir} before
#' launching.
#'
#' @param model_dir Directory containing trained model subdirectories.
#'   Defaults to a "models" folder in the current working directory.
#' @param ... Additional arguments passed to \code{shiny::runApp()}
#' @return Invisible NULL (called for side effect)
#' @export
run_autoSpectra <- function(model_dir = "models", ...) {
  app_dir <- system.file("shiny", package = "autoSpectra")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    # Fallback: run from inst/shiny within the source tree (dev mode)
    app_dir <- file.path(find.package("autoSpectra"), "inst", "shiny")
    if (!dir.exists(app_dir))
      stop("Shiny app directory not found. ",
           "If running from source, ensure inst/shiny/app.R exists.")
  }
  # Pass model_dir as a global option so app.R can read it
  old <- options(autoSpectra.model_dir = model_dir)
  on.exit(options(old))
  shiny::runApp(app_dir, ...)
}
