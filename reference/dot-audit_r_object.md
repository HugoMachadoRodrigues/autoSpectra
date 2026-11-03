# Walk an arbitrary R object and report anything sized like training data

Walk an arbitrary R object and report anything sized like training data

## Usage

``` r
.audit_r_object(
  x,
  path = "",
  max_dim = .AUDIT_MAX_DIM,
  max_chars = .AUDIT_MAX_CHARS,
  max_lines = .AUDIT_MAX_LINES,
  .depth = 0L
)
```

## Arguments

- x:

  Object to inspect

- path:

  Character breadcrumb used in findings

- max_dim:

  Largest permitted extent along any axis

- max_chars:

  Largest permitted length of a single character element

## Value

Character vector of findings (empty when clean)
