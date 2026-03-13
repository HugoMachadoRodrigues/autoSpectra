# Apply Savitzky-Golay filter row-wise to a spectral matrix

Apply Savitzky-Golay filter row-wise to a spectral matrix

## Usage

``` r
apply_sg_matrix(M, m = 11, p = 2, d = 0)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths)

- m:

  Half-window size (total window = 2m+1)

- p:

  Polynomial order

- d:

  Derivative order (0 = smooth, 1 = first derivative)

## Value

Filtered matrix; on error, returns input unchanged
