# Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles, latent stats)

Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles, latent
stats)

## Usage

``` r
save_metrics_json(path_base, metrics, conf_q, feat_range, latent)
```

## Arguments

- path_base:

  File path base (no extension)

- metrics:

  List from metrics_from_y()

- conf_q:

  List with q90 and q95 conformal quantiles

- feat_range:

  Data frame from feature_minmax()

- latent:

  List with mu, Sigma, df, thr95
