# Load a saved soilVAE model and its scaler

Load a saved soilVAE model and its scaler

## Usage

``` r
load_soilVAE(family_id, prop, model_dir = "models")
```

## Arguments

- family_id:

  Family identifier

- prop:

  Soil property name

- model_dir:

  Root model directory (default "models")

## Value

List with `model` (keras) and `scaler` (list with mean, sd)
