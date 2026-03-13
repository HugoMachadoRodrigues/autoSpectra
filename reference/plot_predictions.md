# Plot predicted vs observed values for a soil property

Plot predicted vs observed values for a soil property

## Usage

``` r
plot_predictions(
  observed,
  predicted,
  prop_label = "Property",
  show_metrics = TRUE
)
```

## Arguments

- observed:

  Numeric vector of observed values

- predicted:

  Numeric vector of predicted values

- prop_label:

  Display label for the property (e.g., "Organic Carbon (%)")

- show_metrics:

  Logical; add RMSE, R2, RPIQ annotation (default TRUE)

## Value

A ggplot2 object
