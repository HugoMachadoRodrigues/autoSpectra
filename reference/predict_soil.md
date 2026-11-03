# Predict soil properties from a spectral data frame

Resamples user spectra to the family's canonical grid, applies the
family's preprocessing pipeline (Safanelli's SNV for Cubist; SG-smooth +
masked SNV for masked VAE), generates a per-sample coverage mask from
the resampled `NA` pattern, and dispatches to the appropriate model
loader.

## Usage

``` r
predict_soil(
  df,
  family_id,
  properties = NULL,
  model_dir = getOption("autoSpectra.model_dir", "models"),
  disable_pp = FALSE,
  with_pi = FALSE,
  confidence = 0.68
)
```

## Arguments

- df:

  Data frame with a `Soil_ID` column and numeric spectral columns named
  by wavelength (nm) or wavenumber (cm-1)

- family_id:

  Family identifier from `model_registry`

- properties:

  Properties to predict (default: all)

- model_dir:

  Root directory where models are stored

- disable_pp:

  Skip preprocessing (debug only)

- with_pi:

  Include prediction intervals when available (Cubist family)

- confidence:

  Confidence level for prediction intervals (default 0.68)

## Value

Data frame with `Soil_ID`, one column per predicted property, and (when
`with_pi=TRUE`) `*_lo`/`*_hi` interval columns
