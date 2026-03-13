# Extract OSSL VisNIR spectra as a matrix aligned to the OSSL standard grid

Returns a numeric matrix with rows = samples, columns = wavelengths in
nm (350 to 2500, step 2). The row names are the layer UUIDs.

## Usage

``` r
ossl_visnir_matrix(visnir_df)
```

## Arguments

- visnir_df:

  Raw VisNIR data frame from ossl_load_raw("visnir")

## Value

Numeric matrix (n_samples × 1076 wavelengths)
