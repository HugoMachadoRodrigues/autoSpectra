# Check that a gzipped OSSL cache file is complete

An interrupted download leaves a truncated `.csv.gz` behind. Because the
old cache check was a bare
[`file.exists()`](https://rdrr.io/r/base/files.html), that stub would
then be treated as a valid cache forever, and `fread()` would either
fail or – worse – silently read a partial database.

## Usage

``` r
.ossl_gz_ok(path)
```

## Arguments

- path:

  Path to a `.csv.gz` file

## Value

`TRUE` if the file exists and its gzip stream is complete

## Details

Note that R's own [`gzfile()`](https://rdrr.io/r/base/connections.html)
is NOT sufficient here: it reads a truncated member up to the cut and
then reports a clean EOF, with no error or warning. We therefore verify
the CRC32/ISIZE trailer with `gzip -t`, which is a real gzip
implementation. Where `gzip` is unavailable (some Windows setups) we
fall back to a size sanity check, and rely on the atomic `.part` -\>
rename in
[`ossl_download_file()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_download_file.md)
for correctness.
