# Train soilVAE models for all properties of one model family using OSSL data

This is the main training entry point. It downloads OSSL data (if
needed), applies the family's preprocessing pipeline, and trains one
model per soil property, saving outputs to
`out_dir/<family_id>/models/`.

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

  Family ID from model_registry, e.g. "OSSL_VisNIR"

- out_dir:

  Root output directory (default: "models")

- cache_dir:

  OSSL data cache directory

- properties:

  Properties to train (default: all in family)

- download_if_missing:

  Download OSSL data if not cached

- ...:

  Additional arguments passed to train_soilVAE()

## Value

Invisible named list of results (TRUE/NULL per property)
