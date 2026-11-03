# Canonical Safanelli preprocessing for VisNIR (reflectance input)

Per Safanelli 2025: SNV alone is sufficient and superior to derivatives.
The `ABSORBANCE` step is added for VisNIR because OSSL VisNIR data is
stored as reflectance.

## Usage

``` r
preprocess_visnir(M)
```

## Arguments

- M:

  Numeric matrix of reflectance values

## Value

SNV-preprocessed absorbance matrix
