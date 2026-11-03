# Fit a PCA on a spectral matrix and keep the first n_components scores

Uses
[`irlba::prcomp_irlba()`](https://rdrr.io/pkg/irlba/man/prcomp_irlba.html)
when available (randomised truncated SVD — much less memory and faster
for rank-120 on large matrices like 130k × 1076) and falls back to
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) otherwise.
Returns the fitted object plus centred/scaled metadata for projecting
new spectra and computing Q-statistic residuals.

## Usage

``` r
.fit_pca(X, n_components = 120L)
```

## Arguments

- X:

  Numeric matrix (rows = samples, cols = preprocessed bands)

- n_components:

  Number of PCs to retain (default 120, per Safanelli)

## Value

List with `rotation`, `center`, `scale`, `n_components`, `scores`,
`var_explained`
