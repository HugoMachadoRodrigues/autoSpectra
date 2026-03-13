# Extract OSSL MIR spectra as a matrix aligned to the OSSL standard grid

Returns a numeric matrix with rows = samples, columns = wavenumbers in
cm-1 (600 to 4000, step 2). Values are already in absorbance.

## Usage

``` r
ossl_mir_matrix(mir_df)
```

## Arguments

- mir_df:

  Raw MIR data frame from ossl_load_raw("mir")

## Value

Numeric matrix (n_samples × 1701 wavenumbers)
