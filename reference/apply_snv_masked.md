# Apply SNV using only the masked (covered) bands of each row

For each row, the mean and SD are computed using only the columns where
`mask == 1`. The same row is then z-normalised across all columns
(covered + uncovered), and uncovered positions are zero-filled so the
VAE encoder sees a numerically stable input.

## Usage

``` r
apply_snv_masked(M, mask)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)

- mask:

  Binary matrix of the same shape (1 = covered, 0 = missing)

## Value

SNV matrix; uncovered positions set to 0

## Details

This is the SNV variant used by the masked-VAE families: it makes the
normalisation invariant to which instrument covered which bands.
