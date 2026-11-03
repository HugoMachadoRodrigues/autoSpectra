# Squared Mahalanobis distance against a multivariate Gaussian (mu, Sigma)

Squared Mahalanobis distance against a multivariate Gaussian (mu, Sigma)

## Usage

``` r
mahalanobis_d2(Z, mu, Sigma)
```

## Arguments

- Z:

  Numeric matrix of latent vectors (rows = samples, cols = latent dim)

- mu:

  Numeric vector (length = latent dim)

- Sigma:

  Numeric matrix (latent dim x latent dim)

## Value

Numeric vector of D^2 values, length nrow(Z)
