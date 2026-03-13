# Join OSSL spectra with soil lab data for model training

Join OSSL spectra with soil lab data for model training

## Usage

``` r
ossl_join(spectra_mat, soillab_df, properties = ossl_l1_properties)
```

## Arguments

- spectra_mat:

  Matrix from ossl_visnir_matrix() or ossl_mir_matrix()

- soillab_df:

  Data frame from ossl_soillab()

- properties:

  Properties to include; default all available

## Value

Data frame with Soil_ID + spectral columns + property columns
