# Train a normalised-conformal calibration for a Cubist + PCA model

Mirrors the Safanelli 2025 procedure (eqs. 9-10):

1.  Cross-validated predictions \\\hat{y}\_i^{CV}\\ on the training set
    via internal k-fold CV (Safanelli: "Cross-validated predictions were
    also used to estimate the unbiased absolute error", p.10)

2.  Compute absolute residuals \\\delta_i = \|y_i - \hat{y}\_i^{CV}\|\\

3.  Train an error model with the SAME architecture (Cubist on the same
    PC scores) that maps PC scores to \\\delta_i\\

4.  Compute non-conformity scores \\\alpha_i = \delta_i /
    \hat{\delta}\_i\\

5.  Store the `confidence`-th quantile of \\\alpha_i\\

## Usage

``` r
train_normalized_conformal(
  pca,
  response_model,
  X_train,
  y_train,
  confidence = 0.68,
  committees = NULL,
  k_cv = 5L,
  conformal_subset_n = 10000L,
  seed = 42L
)
```

## Arguments

- pca:

  Fitted PCA list (from
  [`.fit_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-fit_pca.md))

- response_model:

  Cubist model trained on `pca$scores`

- X_train:

  Preprocessed TRAINING spectra used to fit `response_model`

- y_train:

  Numeric vector of training targets

- confidence:

  Confidence level (default 0.68 = ~1 SD per Safanelli)

- committees:

  Cubist committees for the error model (default = same as response
  model when retrievable, else 5)

- k_cv:

  Folds used to obtain CV residuals on the training set (default 5)

- conformal_subset_n:

  Cap on how many training samples feed the inner CV that produces the
  residuals (default 10000). Calibration only needs a stable quantile,
  so a subsample avoids refitting the full training set k_cv times.

- seed:

  RNG seed

## Value

List with `error_model`, `alpha_q`, `confidence`
