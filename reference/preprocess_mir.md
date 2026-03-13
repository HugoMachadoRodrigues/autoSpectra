# Canonical two-step SG preprocessing for MIR (absorbance input)

Canonical two-step SG preprocessing for MIR (absorbance input)

## Usage

``` r
preprocess_mir(M, m = 11, p = 2)
```

## Arguments

- M:

  Numeric matrix of absorbance values (already in absorbance units)

- m:

  Half-window for SG (default 11)

- p:

  Polynomial order (default 2)

## Value

Preprocessed matrix: smooth → 1st derivative
