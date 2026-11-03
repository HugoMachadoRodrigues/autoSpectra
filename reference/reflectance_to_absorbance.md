# Convert reflectance to absorbance

Accepts reflectance in 0-1 range or as percent 0-100 (auto-rescaled).

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
