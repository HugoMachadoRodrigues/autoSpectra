# Download one OSSL component file to the cache directory

Download one OSSL component file to the cache directory

## Usage

``` r
ossl_download_file(
  component = c("visnir", "mir", "soillab"),
  cache_dir = ossl_cache_dir(),
  force = FALSE
)
```

## Arguments

- component:

  One of "visnir", "mir", or "soillab"

- cache_dir:

  Local directory for caching (default: ossl_cache_dir())

- force:

  Re-download even if file already exists

## Value

Invisible path to the downloaded file
