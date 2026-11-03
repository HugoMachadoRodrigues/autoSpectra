# Train soilVAE models for all properties of one model family using OSSL data

Dispatches to `train_soilVAE` or `train_maskedVAE` depending on the
family's `algorithm` field.

## Usage

``` r
train_ossl_models(
  family_id,
  out_dir = "models",
  cache_dir = ossl_cache_dir(),
  properties = NULL,
  download_if_missing = TRUE,
  ...
)
```

## Arguments

- family_id:

  Family ID from `model_registry`

- out_dir:

  Output root (default "models")

- cache_dir:

  OSSL cache directory

- properties:

  Properties to train (default: all in family)

- download_if_missing:

  Auto-download OSSL data if missing

- ...:

  Extra args passed to the per-property trainer

## Value

Invisible named list of per-property results
