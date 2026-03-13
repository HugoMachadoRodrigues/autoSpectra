# Extract OSSL Level-1 soil lab data for a set of properties

Selects the best representative column for each requested property by
picking the column with the fewest NAs among candidates matching the
simplified property name prefix.

## Usage

``` r
ossl_soillab(soillab_df, properties = ossl_l1_properties)
```

## Arguments

- soillab_df:

  Raw soillab data frame from ossl_load_raw("soillab")

- properties:

  Character vector of simplified property names (from
  ossl_l1_properties). Default: all.

## Value

Data frame with columns: id.layer_uuid_txt + one column per property
