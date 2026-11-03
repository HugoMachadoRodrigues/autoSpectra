# Training-scale bounds used to keep a back-transform from extrapolating

Recorded at fit time and stored beside the model, so prediction applies
the same bound the model was calibrated within.

## Usage

``` r
target_bounds(y_model, pad = 0.05)
```

## Arguments

- y_model:

  Numeric vector of TRAINING targets on the modelling scale (already
  log1p-transformed where the policy applies)

- pad:

  Fractional slack on the observed range, so a sample marginally beyond
  the training extremes is not flattened (default 0.05 = 5% of range)

## Value

List with `lo` and `hi` on the modelling scale
