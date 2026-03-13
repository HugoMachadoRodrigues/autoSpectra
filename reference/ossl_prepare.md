# Convenience: load and join all OSSL data for a given sensor type

Downloads data if not yet cached, then returns a joined data frame ready
for model training.

## Usage

``` r
ossl_prepare(
  sensor_type = c("visnir", "mir"),
  cache_dir = ossl_cache_dir(),
  properties = ossl_l1_properties,
  download_if_missing = TRUE
)
```

## Arguments

- sensor_type:

  "visnir" or "mir"

- cache_dir:

  Cache directory

- properties:

  Soil properties to include

- download_if_missing:

  Automatically download if files are absent

## Value

Data frame with Soil_ID + spectral + soil property columns
