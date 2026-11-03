# Compute Q-statistic and in-domain flag for new spectra against a fitted PCA

Compute Q-statistic and in-domain flag for new spectra against a fitted
PCA

## Usage

``` r
predict_q_statistic(pca, X_new, ucl_q)
```

## Arguments

- pca:

  Fitted PCA list (from
  [`.fit_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-fit_pca.md))

- X_new:

  Numeric matrix of new (preprocessed) spectra

- ucl_q:

  Upper Confidence Limit (from
  [`q_statistic_ucl()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/q_statistic_ucl.md))

## Value

Data frame with columns: q_stat, ucl_q, in_domain
