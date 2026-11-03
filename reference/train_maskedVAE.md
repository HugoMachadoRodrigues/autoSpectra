# Train a masked soilVAE for one soil property

Augments the training set via
[`spectral_dropout()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/spectral_dropout.md)
so the encoder learns to handle partial-coverage inputs typical of real
instruments (NeoSpectra, Bruker MPA, PerkinElmer Spectrum One, etc.).

## Usage

``` r
train_maskedVAE(
  X,
  y,
  family_id,
  prop,
  wavegrid,
  sensor_type,
  out_dir = "models",
  epochs = 100L,
  batch_size = 64L,
  val_split = 0.15,
  patience_es = 20L,
  patience_rl = 7L,
  latent_dim = 32L,
  kl_beta = 5e-04,
  kl_warmup_epochs = 30L,
  n_aug_copies = 3L,
  min_n = 50L
)
```

## Arguments

- X:

  Numeric matrix (rows = samples, cols = bands; full-coverage, already
  SNV-normalised on the full grid)

- y:

  Numeric vector of target values

- family_id:

  Family ID (used for output directory)

- prop:

  Soil property name

- wavegrid:

  Numeric vector of wavelength positions (length == ncol(X))

- sensor_type:

  "visnir" or "mir" (controls augmentation cut-offs)

- out_dir:

  Output root (default "models")

- epochs:

  Max epochs (default 80)

- batch_size:

  Mini-batch (default 64)

- val_split:

  Validation fraction (default 0.15)

- patience_es:

  Early-stopping patience (default 10)

- patience_rl:

  LR-reduction patience (default 5)

- latent_dim:

  Latent dim (default 16)

- kl_beta:

  Target KL weight (default 1e-3)

- kl_warmup_epochs:

  Annealing epochs (default 30)

- n_aug_copies:

  Augmented copies per sample (default 2)

- min_n:

  Minimum samples (default 50)

## Value

Invisible TRUE on success, NULL on skip
