# Build Leave-One-Dataset-Out (LODO) k-fold indices

Each fold withholds all samples from one OSSL contributing dataset.
Datasets are randomly assigned to folds when there are more datasets
than folds.

## Usage

``` r
stratified_kfold(dataset_col, k = 10L, seed = 42L)
```

## Arguments

- dataset_col:

  Character vector of dataset labels (one per sample)

- k:

  Number of folds (default 10)

- seed:

  RNG seed

## Value

List of length k of integer test-index vectors
