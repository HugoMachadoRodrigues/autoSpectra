# Check whether pre-trained models are available locally

Check whether pre-trained models are available locally

## Usage

``` r
models_available(
  family_id = c("OSSL_VisNIR", "OSSL_MIR"),
  model_dir = getOption("autoSpectra.model_dir", "models")
)
```

## Arguments

- family_id:

  Family ID (`"OSSL_VisNIR"` or `"OSSL_MIR"`)

- model_dir:

  Root model directory

## Value

Named logical vector: `TRUE` if \>= 5 .h5 files found
