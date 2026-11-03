# Audit a built release archive for redistributed OSSL training data

Extracts the zip and inspects every artefact. Fails on: a file outside
the publishing allowlist, an array or vector with a per-sample axis, or
a character blob large enough to be a serialised training matrix.

## Usage

``` r
audit_release_archive(zip_path, max_dim = .AUDIT_MAX_DIM, quiet = FALSE)
```

## Arguments

- zip_path:

  Path to a `*_models.zip` built by `upload_models.R`

- max_dim:

  Largest permitted extent along any axis (default 5000)

- quiet:

  Suppress the printed report

## Value

Invisibly, a list with `ok`, `files` and `findings`
