# Compute prediction metrics (RMSE, Bias, R², RPIQ, CCC)

Compute prediction metrics (RMSE, Bias, R², RPIQ, CCC)

## Usage

``` r
metrics_from_y(y_true, y_pred)
```

## Arguments

- y_true:

  Numeric vector of observed values

- y_pred:

  Numeric vector of predicted values

## Value

Named list: RMSE, Bias, R2, RPIQ, CCC (Lin's concordance)
