# Train a soilVAE model for a single soil property

Performs train / calibration / test split, fits the model with early
stopping and learning-rate reduction, computes test metrics and
conformal quantiles, and saves the model + scaler + metrics to disk.

## Usage

``` r
train_soilVAE(
  X,
  y,
  family_id,
  prop,
  out_dir = "models",
  epochs = 30L,
  batch_size = 32L,
  val_split = 0.2,
  patience_es = 8L,
  patience_rl = 4L,
  latent_dim = 16L,
  min_n = 30L
)
```

## Arguments

- X:

  Numeric matrix (rows = samples, cols = preprocessed bands)

- y:

  Numeric vector of soil property values (length == nrow(X))

- family_id:

  Family identifier string (used for output directory)

- prop:

  Soil property name (used for output file names)

- out_dir:

  Root output directory (default: "models")

- epochs:

  Max training epochs (default 30)

- batch_size:

  Mini-batch size (default 32)

- val_split:

  Validation fraction from training data (default 0.2)

- patience_es:

  Early stopping patience (default 8)

- patience_rl:

  LR reduction patience (default 4)

- latent_dim:

  Latent dimension (default 16)

- min_n:

  Minimum samples required to attempt training (default 30)

## Value

Invisible TRUE on success, NULL if skipped
