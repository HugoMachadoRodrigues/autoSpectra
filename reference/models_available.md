# Check whether pre-trained models are available locally

Counts both `.keras` (VAE), `.h5` (legacy VAE), and `_cubist.rds`
(Cubist baseline) artefacts.

## Usage

``` r
models_available(
  family_id = c("OSSL_VisNIR", "OSSL_MIR", "OSSL_VisNIR_Cubist", "OSSL_MIR_Cubist"),
  model_dir = getOption("autoSpectra.model_dir", "models")
)
```

## Arguments

- family_id:

  Family ID

- model_dir:

  Root model directory

## Value

Named logical vector: `TRUE` if \>= 5 model artefacts found
