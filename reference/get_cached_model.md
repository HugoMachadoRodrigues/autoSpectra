# Retrieve a soilVAE model from the in-memory cache

On the first call for a given `family_id`/`prop` combination the model
is loaded from disk via
[`load_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/load_soilVAE.md)
and stored in the cache. Subsequent calls return the cached object
instantly.

## Usage

``` r
get_cached_model(
  family_id,
  prop,
  model_dir = getOption("autoSpectra.model_dir", "models")
)
```

## Arguments

- family_id:

  Family identifier (`"OSSL_VisNIR"` or `"OSSL_MIR"`)

- prop:

  Soil property name (e.g. `"oc"`)

- model_dir:

  Root directory where trained models are stored

## Value

List with elements `model` (keras Model) and `scaler` (list with `mean`
and `sd`)
