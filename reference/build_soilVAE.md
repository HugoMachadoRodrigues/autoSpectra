# Build a standard soilVAE

Build a standard soilVAE

## Usage

``` r
build_soilVAE(
  d_in,
  latent_dim = 32L,
  loss_weights = c(0.1, 1),
  l2_reg = 1e-05,
  dropout_rate = 0.2,
  encoder_dropout = 0.05,
  learning_rate = 0.001
)
```

## Arguments

- d_in:

  Integer; number of input spectral bands

- latent_dim:

  Integer; latent space dimension (default 16)

- loss_weights:

  Numeric vector length 2 (rec, pred); default c(0.3, 0.3)

- l2_reg:

  L2 regularisation on Dense layers (default 1e-5)

- dropout_rate:

  Dropout rate for the prediction head (default 0.2)

- encoder_dropout:

  Dropout inside encoder (default 0.05)

- learning_rate:

  Initial Adam learning rate (default 1e-3)

## Value

Compiled keras Model object
