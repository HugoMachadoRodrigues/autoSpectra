# Pre-load all soilVAE models for a family into memory

Iterates over all (or selected) properties for `family_id` and loads
each model into the in-memory cache. Call this once at Shiny app startup
to eliminate per-prediction disk I/O.

## Usage

``` r
preload_ossl_models(
  family_id,
  model_dir = getOption("autoSpectra.model_dir", "models"),
  properties = NULL,
  verbose = TRUE
)
```

## Arguments

- family_id:

  Family identifier (`"OSSL_VisNIR"` or `"OSSL_MIR"`)

- model_dir:

  Root directory where trained models are stored

- properties:

  Character vector of properties to pre-load. `NULL` (default) loads all
  properties defined for the family.

- verbose:

  Print progress messages (default `TRUE`)

## Value

Invisible character vector of successfully loaded property names
