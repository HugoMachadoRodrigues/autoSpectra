# Forward log1p transform for properties whose policy asks for it

Forward log1p transform for properties whose policy asks for it

## Usage

``` r
apply_target_transform(y, prop)
```

## Arguments

- y:

  Numeric vector of raw values

- prop:

  Property name; the transform is applied only when it appears in
  `ossl_l1_log_transform`

## Value

Numeric vector on the modelling scale, with values below -1 set NA
