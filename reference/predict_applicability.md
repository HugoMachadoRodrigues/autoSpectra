# Compute applicability domain score for new samples

Uses the Mahalanobis distance in latent space, calibrated against the
chi-squared threshold at 95% (saved during training).

## Usage

``` r
predict_applicability(df, family_id, prop, model_dir = "models")
```

## Arguments

- df:

  Data frame with Soil_ID and spectral columns

- family_id:

  Family identifier

- prop:

  Soil property name (model must have been trained with latent stats)

- model_dir:

  Root model directory

## Value

Data frame with Soil_ID, mahal_dist, thr95, in_domain (logical)
