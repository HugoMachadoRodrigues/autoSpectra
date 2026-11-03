# Build random 10-fold indices (Safanelli default)

Build random 10-fold indices (Safanelli default)

## Usage

``` r
random_kfold(n, k = 10L, seed = 42L)
```

## Arguments

- n:

  Number of samples

- k:

  Number of folds (default 10)

- seed:

  RNG seed

## Value

List of length k; each element is integer vector of test indices
