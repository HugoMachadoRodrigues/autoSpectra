# Apply a spectral preprocessing pipeline

Step strings:

- `"ABSORBANCE"`: reflectance -\> absorbance (-ln R)

- `"SG_SMOOTH(m,p)"`: Savitzky-Golay smoothing (d=0)

- `"SG_DERIV(m,p,d)"`: Savitzky-Golay derivative (d\>=1)

- `"SNV"`: Standard Normal Variate (Safanelli 2025)

- `"SNV_MASKED"`: per-mask SNV; requires `mask` argument

## Usage

``` r
apply_pipeline(M, steps, mask = NULL, absorbance_base10 = FALSE)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths/wavenumbers)

- steps:

  Character vector of preprocessing step strings

- mask:

  Optional binary matrix for `SNV_MASKED` (same shape as M). If omitted
  while a `SNV_MASKED` step is present, a fully-covered mask is assumed.

- absorbance_base10:

  Logical; use base-10 log for absorbance conversion

## Value

Preprocessed matrix
