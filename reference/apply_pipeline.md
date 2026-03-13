# Apply a spectral preprocessing pipeline

The pipeline is defined as a character vector of step names:

- `"ABSORBANCE"`: convert reflectance → absorbance (−ln R)

- `"SG_SMOOTH(m,p)"`: Savitzky-Golay smooth (derivative = 0)

- `"SG_DERIV(m,p,d)"`: Savitzky-Golay derivative (d ≥ 1)

- `"SG(m,p,d)"`: legacy single-step SG (backward compatible)

- `"SNV"`: Standard Normal Variate

## Usage

``` r
apply_pipeline(M, steps, absorbance_base10 = FALSE)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)

- steps:

  Character vector of preprocessing step strings

- absorbance_base10:

  Logical; use base-10 log for absorbance conversion

## Value

Preprocessed matrix

## Details

The recommended two-step pipeline for OSSL VisNIR data is:
`c("ABSORBANCE", "SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)")`

For MIR data (already in absorbance):
`c("SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)")`
