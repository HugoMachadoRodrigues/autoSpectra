# Invert the target transform without letting it explode

Clamps to the calibration range on the modelling scale BEFORE
exponentiating, then inverts. For untransformed properties the clamp
still applies, which costs nothing and keeps predictions physically
plausible.

## Usage

``` r
invert_target_transform(y_hat, prop, bounds = NULL, floor_zero = TRUE)
```

## Arguments

- y_hat:

  Numeric vector of predictions on the modelling scale

- prop:

  Property name

- bounds:

  List from
  [`target_bounds()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/target_bounds.md),
  or `NULL` to skip clamping (which restores the old, unbounded
  behaviour)

- floor_zero:

  Clamp the result at 0 for non-negative properties (default `TRUE`;
  concentrations cannot be negative)

## Value

Numeric vector in original units
