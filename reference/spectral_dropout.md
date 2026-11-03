# Random spectral-dropout augmentation for masked-VAE training

For every input row, generates additional copies that simulate partial-
coverage instruments by zeroing out the mask outside randomised cut-off
wavelengths. Per OSSL v1.2, the realistic VisNIR cut-offs are 350, 400,
800, 1000, 1350 nm (corresponding to ASD, Foss, Bruker MPA, PerkinElmer,
NeoSpectra); MIR cut-offs cluster around 600-700 cm-1.

## Usage

``` r
spectral_dropout(
  X,
  wavegrid,
  sensor_type = c("visnir", "mir"),
  n_copies = 2L,
  seed = NULL
)
```

## Arguments

- X:

  Numeric matrix of preprocessed spectra (rows = samples, cols = bands)

- wavegrid:

  Numeric vector of wavelength/wavenumber positions (length == ncol(X))

- sensor_type:

  "visnir" or "mir" (controls realistic cut-off pool)

- n_copies:

  How many augmented copies to generate per input row (in addition to
  the original full-coverage copy). Default 2.

- seed:

  Optional RNG seed for reproducibility

## Value

List with `X` (rbind of all copies) and `mask` (same shape, 0/1)
