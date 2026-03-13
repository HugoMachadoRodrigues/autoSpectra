# Plot a mean spectrum with ±1 SD ribbon from a spectral matrix

Plot a mean spectrum with ±1 SD ribbon from a spectral matrix

## Usage

``` r
plot_mean_spectrum(
  M,
  wl,
  xlab = "Wavelength (nm)",
  ylab = "Mean Response",
  title = "Mean ± SD spectrum",
  colour = "#2166ac"
)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = wavelengths)

- wl:

  Numeric vector of wavelength positions (same length as ncol(M))

- xlab:

  X-axis label (default "Wavelength (nm)")

- ylab:

  Y-axis label (default "Mean Response")

- title:

  Plot title

- colour:

  Line/ribbon colour (default "#2166ac")

## Value

A ggplot2 object
