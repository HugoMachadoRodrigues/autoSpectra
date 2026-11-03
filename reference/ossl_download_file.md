# Download one OSSL component file to the cache directory

Downloads from the official OSSL distribution bucket documented at
<https://docs.soilspectroscopy.org/db-access.html>. autoSpectra never
bundles or re-hosts OSSL data — each user fetches it from the official
source into their own local cache.

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

  Re-download even if a valid file already exists

## Value

Invisible path to the downloaded file

## Details

The file is written to a temporary path and only moved into place once
the HTTP status and the gzip stream both check out, so an interrupted
download can never masquerade as a valid cache entry.
