# Launch the autoSpectra Shiny application

Opens the interactive soil spectral prediction interface in the default
browser. Models must have been trained and saved to `model_dir` before
launching.

## Usage

``` r
run_autoSpectra(model_dir = "models", ...)
```

## Arguments

- model_dir:

  Directory containing trained model subdirectories. Defaults to a
  "models" folder in the current working directory.

- ...:

  Additional arguments passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html)

## Value

Invisible NULL (called for side effect)
