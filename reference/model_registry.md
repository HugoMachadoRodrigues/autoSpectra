# Official autoSpectra model family registry

Four sensor-agnostic families trained on the full OSSL v1.2 corpus. Two
model types per spectral domain:

- soilVAE (masked, sensor-agnostic) - production models

- Cubist + PCA(120) (Safanelli 2025 replication) - benchmark models

## Usage

``` r
model_registry
```

## Details

All families share the same OSSL spectral grids: VisNIR: 350-2500 nm at
2 nm MIR : 600-4000 cm-1 at 2 cm-1
