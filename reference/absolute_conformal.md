# Simple absolute-residual conformal calibration (used by the VAE family)

Less sophisticated than the normalised variant, but does not require a
second model. Returns a single half-width per confidence level.

## Usage

``` r
absolute_conformal(y_true, y_pred, confidence = c(0.68, 0.9, 0.95))
```

## Arguments

- y_true:

  Calibration observed values

- y_pred:

  Calibration predicted values

- confidence:

  Numeric vector of confidence levels (default c(0.68, 0.90, 0.95))

## Value

Named numeric vector of half-widths (one per confidence level)
