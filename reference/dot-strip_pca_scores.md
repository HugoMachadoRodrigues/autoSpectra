# Drop the per-sample PC scores from a fitted PCA before serialisation

`pca$scores` holds the 120 PC scores of every OSSL training sample and
[`.reconstruct_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-reconstruct_pca.md)
inverts them back to the input spectra at ~99.98% of the variance. Only
the projection basis (rotation/center/scale) is needed at inference time
—
[`.project_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-project_pca.md)
never reads `scores`.

## Usage

``` r
.strip_pca_scores(pca)
```

## Arguments

- pca:

  A list from
  [`.fit_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-fit_pca.md)

## Value

The same list with `scores` replaced by a row count
