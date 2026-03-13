# Compute safe SG parameters given available number of columns

Compute safe SG parameters given available number of columns

## Usage

``` r
sg_safe_params(m, p, ncols)
```

## Arguments

- m:

  Half-window size (prospectr convention)

- p:

  Polynomial order

- ncols:

  Number of spectral columns

## Value

Adjusted list(m, p), or NULL if impossible
