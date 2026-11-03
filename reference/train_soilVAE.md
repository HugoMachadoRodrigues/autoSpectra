# Train a (standard) soilVAE model for one soil property

Train a (standard) soilVAE model for one soil property

## Usage

``` r
train_soilVAE(
  X,
  y,
  family_id,
  prop,
  out_dir = "models",
  epochs = 100L,
  batch_size = 64L,
  val_split = 0.15,
  patience_es = 20L,
  patience_rl = 7L,
  latent_dim = 32L,
  kl_beta = 5e-04,
  kl_warmup_epochs = 30L,
  min_n = 50L
)
```

## Arguments

- X:

  Numeric matrix of preprocessed spectra (full coverage)

- y:

  Numeric vector of target values (length == nrow(X))

- family_id:

  Family ID (e.g. "OSSL_VisNIR_Cubist" - not used for VAE)

- prop:

  Soil property name (used for output filenames)

- out_dir:

  Root output directory (default "models")

- epochs:

  Max epochs (default 80)

- batch_size:

  Mini-batch size (default 32)

- val_split:

  Validation fraction (default 0.2)

- patience_es:

  Early-stopping patience (default 10)

- patience_rl:

  LR-reduction patience (default 5)

- latent_dim:

  Latent dim (default 16)

- kl_beta:

  Target KL weight after annealing (default 1e-3)

- kl_warmup_epochs:

  Annealing warmup (default 30)

- min_n:

  Minimum samples to attempt training (default 50)

## Value

Invisible TRUE on success, NULL if skipped
