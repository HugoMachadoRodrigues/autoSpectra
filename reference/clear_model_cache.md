# Clear the in-memory model cache

Releases all cached keras model objects. Models will be reloaded from
disk on the next prediction call.

## Usage

``` r
clear_model_cache()
```

## Value

Invisible `NULL`
