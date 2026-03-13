# Build the soilVAE asymmetric autoencoder model

Requires the keras package. Call `keras::use_backend("tensorflow")`
before training if needed.

## Usage

``` r
build_soilVAE(d_in, latent_dim = 16L, loss_weights = c(0.3, 0.3))
```

## Arguments

- d_in:

  Integer; number of input features (preprocessed spectral bands)

- latent_dim:

  Integer; latent space dimension (default 16)

- loss_weights:

  Numeric vector length 2; weights for reconstruction and prediction
  losses (default c(0.3, 0.3))

## Value

A compiled keras Model object
