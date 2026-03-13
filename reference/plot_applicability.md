# Plot latent-space applicability domain scores

Mahalanobis distances are shown as a horizontal bar chart; samples
outside the 95% threshold are highlighted in red.

## Usage

``` r
plot_applicability(app_df, title = "Applicability Domain")
```

## Arguments

- app_df:

  Data frame from predict_applicability()

- title:

  Plot title

## Value

A ggplot2 object
