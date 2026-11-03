# Save a z-score scaler to RDS and JSON

Save a z-score scaler to RDS and JSON

## Usage

``` r
save_scaler(mean_y, sd_y, path_base, bounds = NULL)
```

## Arguments

- mean_y:

  Mean of the target variable

- sd_y:

  Standard deviation of the target variable

- path_base:

  File path base (no extension)

- bounds:

  Optional list from
  [`target_bounds()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/target_bounds.md):
  the calibration range on the modelling scale, stored so prediction can
  bound the back-transform the same way training scored it. `NULL`
  stores no range.
