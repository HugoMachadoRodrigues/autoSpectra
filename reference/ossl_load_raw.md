# Load a cached OSSL CSV (gzipped) as a data.table

Load a cached OSSL CSV (gzipped) as a data.table

## Usage

``` r
ossl_load_raw(
  component = c("visnir", "mir", "soillab"),
  cache_dir = ossl_cache_dir()
)
```

## Arguments

- component:

  "visnir", "mir", or "soillab"

- cache_dir:

  Cache directory

## Value

data.table
