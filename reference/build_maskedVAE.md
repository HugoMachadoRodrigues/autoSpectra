# Build a masked soilVAE for sensor-agnostic prediction

Input layout: a single tensor of shape `(2 * d_spec)` where the first
`d_spec` columns are the SNV-normalised spectrum (uncovered bands
zero-filled) and the next `d_spec` columns are the binary coverage mask
(1 = covered, 0 = missing).

## Usage

``` r
build_maskedVAE(
  d_spec,
  latent_dim = 32L,
  loss_weights = c(0.1, 1),
  l2_reg = 1e-05,
  dropout_rate = 0.2,
  encoder_dropout = 0.05,
  learning_rate = 0.001
)
```

## Arguments

- d_spec:

  Integer; number of bands per channel (e.g., 1076 for VisNIR)

- latent_dim:

  Integer; latent dim (default 16)

- loss_weights:

  Numeric vector length 2 (rec, pred); default c(0.3, 0.3)

- l2_reg:

  L2 regularisation (default 1e-5)

- dropout_rate:

  Prediction-head dropout (default 0.2)

- encoder_dropout:

  Encoder dropout (default 0.05)

- learning_rate:

  Initial Adam LR (default 1e-3)

## Value

Compiled keras Model

## Details

Reconstruction loss is masked: only covered bands contribute. The
prediction head reads only from the latent z, so its output is invariant
to which bands were missing.
