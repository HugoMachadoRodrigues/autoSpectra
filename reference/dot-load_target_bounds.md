# Read a model's calibration range, if it was recorded at training time

Models trained before the bounded back-transform have no stored range;
they return NULL and fall back to the old unbounded behaviour rather
than failing.

## Usage

``` r
.load_target_bounds(family_id, prop, model_dir)
```
