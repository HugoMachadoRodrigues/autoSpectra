# Perform 10-fold CV with refitting for one (family, property) combination

Dispatches to either VAE or Cubist+PCA training on each fold based on
the family's algorithm.

## Usage

``` r
cv_train_property(X, y, family_id, prop, folds, fam, ...)
```

## Arguments

- X:

  Numeric matrix of preprocessed (SNV) spectra

- y:

  Numeric target vector

- family_id:

  Family identifier

- prop:

  Soil property name

- folds:

  List of integer test-index vectors (from
  [`random_kfold()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/random_kfold.md)
  or
  [`stratified_kfold()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/stratified_kfold.md))

- fam:

  Family list from `model_registry`

- ...:

  Extra args forwarded to per-fold trainers

## Value

List with per-fold metrics and aggregated mean/SD
