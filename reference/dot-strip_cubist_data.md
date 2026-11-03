# Strip embedded training data from a fitted Cubist model

[`Cubist::cubist()`](http://topepo.github.io/Cubist/reference/cubist.default.md)
keeps a verbatim copy of the whole training matrix (all predictor
columns plus the response) in `$data`, as a CSV text blob. Because we
train on PCA scores of OSSL spectra, that blob would let anyone
reconstruct the OSSL spectra (via
[`.reconstruct_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dot-reconstruct_pca.md))
and read the OSSL lab values straight out of a published model file —
i.e. it would turn our model archive into an unlicensed redistribution
of the OSSL database. It is also pure dead weight: predictions are
bit-identical without it.

## Usage

``` r
.strip_cubist_data(mdl)
```

## Arguments

- mdl:

  A model returned by
  [`Cubist::cubist()`](http://topepo.github.io/Cubist/reference/cubist.default.md)

## Value

The same model with `$data` emptied

## Details

Always call this before serialising a Cubist model.
