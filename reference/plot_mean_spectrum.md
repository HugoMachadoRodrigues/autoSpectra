# Plot a mean spectrum with +/-1 SD ribbon

Accepts either a spectral data frame (with `Soil_ID` + numeric
wavelength columns) plus a `family` from the model registry, OR a raw
matrix plus an explicit wavelength vector.

## Usage

``` r
plot_mean_spectrum(
  x,
  wl = NULL,
  family = NULL,
  id_col = "Soil_ID",
  xlab = NULL,
  ylab = "Mean Response",
  title = "Mean ± SD spectrum",
  colour = "#2166ac"
)
```

## Arguments

- x:

  Either a data frame with `Soil_ID` + numeric spectral columns, OR a
  numeric matrix (rows = samples, cols = wavelengths). The legacy
  argument name `M` is retained as an alias.

- wl:

  Numeric vector of wavelength positions. Required when `x` is a matrix;
  auto-detected when `x` is a data frame.

- family:

  Optional family list from `model_registry`; when provided with a data
  frame, draws the model grid range as dashed vertical lines.

- id_col:

  Sample identifier column name (default `"Soil_ID"`); only used when
  `x` is a data frame.

- xlab:

  X-axis label; auto-detected when `NULL`

- ylab:

  Y-axis label (default `"Mean Response"`)

- title:

  Plot title

- colour:

  Line/ribbon colour (default `"#2166ac"`)

## Value

A ggplot2 object
