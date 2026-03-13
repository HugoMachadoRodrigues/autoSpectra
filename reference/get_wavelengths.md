# Extract wavelength/wavenumber positions from a data frame's column names

Extract wavelength/wavenumber positions from a data frame's column names

## Usage

``` r
get_wavelengths(df, id_col = "Soil_ID")
```

## Arguments

- df:

  A data frame with spectral columns named as numbers (nm or cm-1)

- id_col:

  Name of the sample identifier column to exclude

## Value

A list with `wl` (numeric positions) and `cols` (column names)
