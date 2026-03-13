# Predict soil properties from a spectral data frame using soilVAE models

Predict soil properties from a spectral data frame using soilVAE models

## Usage

``` r
predict_soil(
  df,
  family_id,
  properties = NULL,
  model_dir = getOption("autoSpectra.model_dir", "models"),
  disable_pp = FALSE
)
```

## Arguments

- df:

  Data frame with a `Soil_ID` column and numeric spectral columns named
  by wavelength (nm) or wavenumber (cm-1)

- family_id:

  Family identifier: `"OSSL_VisNIR"` or `"OSSL_MIR"`

- properties:

  Character vector of OSSL L1 property keys to predict. Default (`NULL`)
  predicts all properties available for the family.

- model_dir:

  Root directory where trained models are stored

- disable_pp:

  Logical; skip spectral preprocessing (debug only)

## Value

Data frame with `Soil_ID` and one column per predicted property
