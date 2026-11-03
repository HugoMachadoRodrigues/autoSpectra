# Is this artefact cleared for publication?

Applied by `upload_models.R` to every file before it is zipped for
Zenodo. Anything unmatched aborts the upload so a new artefact must be
reviewed deliberately rather than shipped by default.

## Usage

``` r
.is_publishable(paths)
```

## Arguments

- paths:

  Character vector of file paths (only the basename is tested)

## Value

Logical vector, `TRUE` where the file may be published
