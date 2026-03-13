# Convert reflectance to absorbance

Accepts reflectance in the range 0-1 or as percentage 0-100. Values in
the 0-100 range are automatically rescaled.

## Usage

``` r
reflectance_to_absorbance(M, base10 = FALSE)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths)

- base10:

  Logical; if TRUE use -log10(R), otherwise -ln(R)

## Value

Absorbance matrix with same dimensions
