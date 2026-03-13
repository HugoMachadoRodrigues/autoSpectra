# Resample a spectral matrix to a target wavelength grid via linear interpolation

Resample a spectral matrix to a target wavelength grid via linear
interpolation

## Usage

``` r
resample_to_grid(M, src_wl, target_wl)
```

## Arguments

- M:

  Numeric matrix (rows = samples, cols = source wavelengths)

- src_wl:

  Numeric vector of source wavelength positions

- target_wl:

  Numeric vector of target wavelength positions

## Value

Numeric matrix resampled to target_wl
