# Project new spectra into a fitted PCA space

Project new spectra into a fitted PCA space

## Usage

``` r
.project_pca(pca, X_new)
```

## Arguments

- pca:

  Object returned by
  [`.fit_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-fit_pca.md)

- X_new:

  Numeric matrix in the same band order as the fitting set

## Value

Matrix of PC scores (n_new x n_components)
