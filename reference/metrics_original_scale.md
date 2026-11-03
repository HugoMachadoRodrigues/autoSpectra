# Score predictions in original units, with the bounded back-transform

The single place training code should compute metrics for a transformed
property, so every family scores the same way.

## Usage

``` r
metrics_original_scale(y_true_model, y_pred_model, prop, bounds = NULL)
```

## Arguments

- y_true_model:

  Observed values on the modelling scale

- y_pred_model:

  Predicted values on the modelling scale

- prop:

  Property name

- bounds:

  List from
  [`target_bounds()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/target_bounds.md)
  (typically built from the training targets), or `NULL`

## Value

List of metrics from
[`metrics_from_y()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/metrics_from_y.md),
in original units
