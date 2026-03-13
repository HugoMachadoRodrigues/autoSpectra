# Plot one or more spectra from a data frame

Plot one or more spectra from a data frame

## Usage

``` r
plot_spectra(
  df,
  sample_ids = NULL,
  id_col = "Soil_ID",
  xlab = NULL,
  ylab = "Response",
  colour_by = NULL,
  family = NULL,
  title = "Spectra",
  alpha = 0.7
)
```

## Arguments

- df:

  Data frame with Soil_ID column and numeric spectral columns

- sample_ids:

  Character vector of Soil_IDs to plot. If NULL, plots all.

- id_col:

  Name of the sample ID column (default "Soil_ID")

- xlab:

  X-axis label; auto-detected if NULL ("Wavelength (nm)" or "Wavenumber
  (cm-1)")

- ylab:

  Y-axis label (default "Response")

- colour_by:

  Name of a non-spectral column to colour lines by; NULL for a single
  colour

- family:

  Family list (from model_registry) used to draw the model grid range as
  dashed vertical lines. Pass NULL to skip.

- title:

  Plot title (default: "Spectra")

- alpha:

  Line transparency (default 0.7)

## Value

A ggplot2 object
