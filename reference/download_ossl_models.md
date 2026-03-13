# Download pre-trained OSSL soilVAE models from Zenodo

Downloads one or both pre-trained model archives from the autoSpectra
Zenodo deposit and extracts them into `model_dir`. After the first
download, models are cached locally and loaded from disk (or from the
in-memory cache via
[`get_cached_model()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/get_cached_model.md))
on subsequent calls to
[`predict_soil()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_soil.md).

## Usage

``` r
download_ossl_models(
  family_id = c("OSSL_VisNIR", "OSSL_MIR"),
  model_dir = getOption("autoSpectra.model_dir", "models"),
  zenodo_record = .MODELS_ZENODO_RECORD,
  overwrite = FALSE,
  timeout_sec = 3600L
)
```

## Arguments

- family_id:

  Character vector of family IDs to download. One or both of
  `"OSSL_VisNIR"` and `"OSSL_MIR"`. Default: both families.

- model_dir:

  Local directory where models are stored. Default:
  `getOption("autoSpectra.model_dir", "models")`.

- zenodo_record:

  Zenodo record ID (string). Normally left at default.

- overwrite:

  Logical. If `TRUE`, re-download even if models are already present.
  Default `FALSE`.

- timeout_sec:

  Download timeout in seconds (default 3600 = 1 h).

## Value

Invisible character vector of directories written.

## Quick-start workflow

    library(autoSpectra)
    download_ossl_models()                  # one-time setup
    predict_soil(df, "OSSL_VisNIR")         # works immediately
