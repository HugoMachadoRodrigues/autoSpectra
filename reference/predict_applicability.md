# Compute applicability-domain scores for new spectra

Dispatches to Q-statistic (Cubist family) or Mahalanobis on latent space
(VAE family) based on the family's algorithm.

## Usage

``` r
predict_applicability(
  df,
  family_id,
  prop,
  model_dir = getOption("autoSpectra.model_dir", "models")
)
```

## Arguments

- df:

  Data frame with `Soil_ID` + spectral columns

- family_id:

  Family identifier

- prop:

  Soil property name (latent statistics are per-property for VAE)

- model_dir:

  Root model directory

## Value

Data frame with columns: `Soil_ID`, `score`, `threshold`, `in_domain`,
plus `method` attribute
