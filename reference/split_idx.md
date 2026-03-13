# Split indices into train / calibration / test sets

Split indices into train / calibration / test sets

## Usage

``` r
split_idx(n, seed = 42, p_train = 0.6, p_cal = 0.2, p_test = 0.2)
```

## Arguments

- n:

  Total number of samples

- seed:

  Random seed for reproducibility

- p_train:

  Fraction for training

- p_cal:

  Fraction for calibration (conformal prediction)

- p_test:

  Fraction for test

## Value

Named list with integer index vectors: train, calib, test
