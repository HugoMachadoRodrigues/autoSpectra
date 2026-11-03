# Which rows carry any spectral signal at all

The OSSL L0 tables have one row per soil LAYER, not per scan: a layer
with no scan on a given instrument still gets a row, with every spectral
column empty. That is 52.4% of the VisNIR table and 36.9% of MIR.

## Usage

``` r
has_spectral_coverage(mask)
```

## Arguments

- mask:

  Integer/logical coverage matrix from
  [`mask_from_coverage()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/mask_from_coverage.md)

## Value

Logical vector, `TRUE` for rows with at least one covered band

## Details

Such rows MUST be dropped before training. They cannot be caught
downstream with `is.finite(rowSums(X))`, because the pipeline replaces
non-finite cells with 0 first: an all-NA row becomes an all-zero row,
passes [`is.finite()`](https://rdrr.io/r/base/is.finite.html), and
survives into training as a constant vector (constant -1 after SNV).
Every such row is then identical to every other, so identical inputs
force identical predictions and impose a hard, model-independent ceiling
on R2 – which is why two unrelated algorithms previously landed within
0.05 of each other on every property.

Partial coverage is fine, and is precisely what the masked soilVAE
exists to handle; only rows with no coverage at all are unusable.

A minimum-band threshold was considered and rejected on the data:
coverage in OSSL v1.2 is perfectly bimodal. Every non-empty row has at
least 1,051 of 1,076 VisNIR bands (the 1,051 floor is LUCAS/Metrohm
starting at 400 nm rather than 350) or all 1,701 MIR bands. Nothing sits
in between, so `> 0` and `>= 500` select exactly the same rows and a
threshold would be arbitrary complexity. Re-check this if a future OSSL
release mixes in instruments with genuinely sparse coverage.
