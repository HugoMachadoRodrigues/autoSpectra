# Compute applicability-domain scores for new samples

Uses the squared Mahalanobis distance in the soilVAE latent space,
compared against the chi-squared threshold at df = 16, alpha = 0.05
(thr95 ~26.3). Latent statistics (mu, Sigma) are read from the
per-property metrics JSON saved during training.

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

  Data frame with `Soil_ID` and spectral columns

- family_id:

  Family identifier

- prop:

  Soil property name (determines which model's latent space is used)

- model_dir:

  Root model directory

## Value

Data frame with columns `Soil_ID`, `mahal_dist`, `thr95`, `in_domain`
