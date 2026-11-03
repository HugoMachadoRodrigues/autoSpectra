# Plot applicability-domain scores

Works with the unified data frame returned by
[`predict_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_applicability.md)
for both algorithms. Cubist family uses Q-statistic; VAE family uses
Mahalanobis distance in latent space. The plot's x-axis label and
caption adapt to the method recorded in `attr(app_df, "method")`.

## Usage

``` r
plot_applicability(app_df, title = "Applicability Domain")
```

## Arguments

- app_df:

  Data frame from
  [`predict_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_applicability.md)
  with columns `Soil_ID`, `score`, `threshold`, `in_domain`

- title:

  Plot title

## Value

A ggplot2 object
