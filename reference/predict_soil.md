# Predict soil properties from a spectral data frame using saved soilVAE models

Predict soil properties from a spectral data frame using saved soilVAE
models

## Usage

``` r
predict_soil(
  df,
  family_id,
  properties = NULL,
  model_dir = "models",
  disable_pp = FALSE
)
```

## Arguments

- df:

  Data frame with a Soil_ID column and spectral columns (numeric names
  in nm or cm-1)

- family_id:

  Family identifier (e.g., "OSSL_VisNIR", "ASD_DRY")

- properties:

  Character vector of properties to predict. Default: all available in
  the family.

- model_dir:

  Root directory where trained models are stored (default "models")

- disable_pp:

  Logical; skip preprocessing (for debug only)

## Value

Data frame with Soil_ID and one column per predicted property
