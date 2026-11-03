# Train one Cubist + PCA(120) model for a single soil property

Faithful to Safanelli 2025: PCA(120) on SNV-preprocessed spectra,
optional log1p target transform, Cubist with inner-CV committee tuning,
Q-statistic UCL on training residuals, normalised-conformal calibration
for prediction intervals.

## Usage

``` r
train_cubist_pca(
  X,
  y,
  family_id,
  prop,
  out_dir = "models",
  n_components = 120L,
  committees_grid = .CUBIST_COMMITTEE_GRID,
  k_inner = 5L,
  inner_subset_n = 2000L,
  min_n = 50L
)
```

## Arguments

- X:

  Numeric matrix of preprocessed (SNV) spectra

- y:

  Numeric target vector

- family_id:

  Family ID (e.g. "OSSL_VisNIR_Cubist")

- prop:

  Soil property name

- out_dir:

  Output root (default "models")

- n_components:

  PCA components (default 120)

- committees_grid:

  Cubist committee values to tune

- k_inner:

  Inner-CV folds for tuning (default 5)

- inner_subset_n:

  Subsample for tuning (default 2000)

- min_n:

  Minimum samples to attempt training (default 50)

## Value

Invisible TRUE on success, NULL on skip
