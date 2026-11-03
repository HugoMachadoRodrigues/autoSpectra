# Resample a spectral matrix to a target wavelength grid via linear interpolation

Rows with fewer than 2 finite values (e.g., instruments that do not
cover the target range) are left as `NA` in the output and will be
excluded by the `is.finite(rowSums(X))` filter in the training pipeline.

## Usage

``` r
resample_to_grid(M, src_wl, target_wl)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = source wavelengths)

- src_wl:

  Numeric vector of source wavelength positions

- target_wl:

  Numeric vector of target wavelength positions

## Value

Numeric matrix resampled to target_wl
