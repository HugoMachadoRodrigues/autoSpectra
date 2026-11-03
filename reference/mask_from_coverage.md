# Build a coverage mask from a resampled spectral matrix

After
[`resample_to_grid()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/resample_to_grid.md)
samples that the source instrument did not cover end up as `NA`
(`rule = 1` in
[`stats::approx`](https://rdrr.io/r/stats/approxfun.html)). This helper
turns that `NA` pattern into a 0/1 mask.

## Usage

``` r
mask_from_coverage(M_resampled)
```

## Arguments

- M_resampled:

  Numeric matrix returned by
  [`resample_to_grid()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/resample_to_grid.md)

## Value

Integer matrix of the same shape (1 = covered, 0 = missing/NA)
