# Predict with prediction intervals using a normalised-conformal calibration

Predict with prediction intervals using a normalised-conformal
calibration

## Usage

``` r
predict_with_pi(pca, response_model, error_model, alpha_q, X_new)
```

## Arguments

- pca:

  Fitted PCA

- response_model:

  Cubist response model

- error_model:

  Cubist error model from
  [`train_normalized_conformal()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_normalized_conformal.md)

- alpha_q:

  Calibrated non-conformity quantile (alpha_q)

- X_new:

  Preprocessed spectra (matrix)

## Value

Data frame with columns: y_pred, delta_hat, pi_lo, pi_hi
