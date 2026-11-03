# Compute per-band min/max feature range from a preprocessed matrix

Falls back to positional band ids when the matrix has no column names –
otherwise [`colnames()`](https://rdrr.io/r/base/colnames.html) returns
NULL and [`data.frame()`](https://rdrr.io/r/base/data.frame.html) fails
with "arguments imply differing number of rows: 0, 2".

## Usage

``` r
feature_minmax(M)
```

## Arguments

- M:

  Numeric matrix, ideally with named columns

## Value

Data frame with columns: band, min, max
