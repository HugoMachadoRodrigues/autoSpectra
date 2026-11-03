# Load a saved Cubist + PCA bundle for a property

Load a saved Cubist + PCA bundle for a property

## Usage

``` r
load_cubist_pca(family_id, prop, model_dir = "models")
```

## Arguments

- family_id:

  Family ID

- prop:

  Soil property name

- model_dir:

  Root model directory (default "models")

## Value

List with `model` (Cubist), `pca`, `best_committees`
