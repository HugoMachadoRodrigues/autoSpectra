# Inner 5-fold CV to choose the best Cubist committees value

Per Safanelli: hyperparameter search on a SUBSET of ~2000 samples for
speed. We mirror that with `inner_subset_n = 2000`.

## Usage

``` r
.tune_cubist_committees(
  T_scores,
  y,
  committees_grid = .CUBIST_COMMITTEE_GRID,
  k_inner = 5L,
  inner_subset_n = 2000L,
  seed = 42L
)
```

## Arguments

- T_scores:

  Matrix of PC scores (rows = samples, cols = PCs)

- y:

  Numeric target vector

- committees_grid:

  Integer vector of committee values to try

- k_inner:

  Inner-CV folds (default 5, per Safanelli)

- inner_subset_n:

  Subsample size for tuning (default 2000)

- seed:

  RNG seed

## Value

List with `best_committees` and per-config RMSE means
