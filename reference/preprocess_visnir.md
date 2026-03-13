# Canonical two-step SG preprocessing for VisNIR (reflectance input)

Canonical two-step SG preprocessing for VisNIR (reflectance input)

## Usage

``` r
preprocess_visnir(M, m = 11, p = 2)
```

## Arguments

- M:

  Numeric matrix of reflectance values

- m:

  Half-window for SG (default 11)

- p:

  Polynomial order (default 2)

## Value

Preprocessed matrix: absorbance → smooth → 1st derivative
