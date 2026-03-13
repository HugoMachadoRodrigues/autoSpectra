# Download all three OSSL components (VisNIR, MIR, soillab)

Download all three OSSL components (VisNIR, MIR, soillab)

## Usage

``` r
ossl_download(
  cache_dir = ossl_cache_dir(),
  force = FALSE,
  components = c("visnir", "mir", "soillab")
)
```

## Arguments

- cache_dir:

  Local directory for caching

- force:

  Re-download even if files already exist

- components:

  Which components to download (default: all three)

## Value

Invisible named list of file paths
