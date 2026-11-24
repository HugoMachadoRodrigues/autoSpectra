# autoSpectra <img src="man/figures/logo.png" align="right" height="138"/>

<!-- badges: start -->
[![Version](https://img.shields.io/badge/version-0.5.0-blue.svg)](https://github.com/HugoMachadoRodrigues/autoSpectra/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19004686.svg)](https://doi.org/10.5281/zenodo.19004686)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg?logo=R)](https://www.r-project.org/)
[![OSSL](https://img.shields.io/badge/data-OSSL%20v1.2-2e7d32.svg)](https://docs.soilspectroscopy.org)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8070--8126-A6CE39.svg?logo=orcid)](https://orcid.org/0000-0002-8070-8126)
[![R-CMD-check](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/pkgdown.yml/badge.svg)](https://HugoMachadoRodrigues.github.io/autoSpectra/)
<!-- badges: end -->

> **Soil spectral modeling, visualization, and prediction** — models trained on the Open Soil Spectral Library (OSSL v1.2) for VisNIR and MIR spectra: a Cubist + PCA replication of Safanelli et al. (2025), a PCA + soilVAE ensemble that takes exactly the same input and is more accurate, and PCA-VAE models for single instrument windows, so that spectra covering only part of the range are predicted too (`predict_soil_agnostic()` picks the window). `train_sensor_model()` builds the same kind of model for your own instrument.

---

## Table of Contents

1. [Overview](#overview)
2. [Models](#models)
3. [Methodology](#methodology--what-we-copied-from-safanelli-2025-and-what-we-added)
4. [Supported Instruments](#supported-instruments)
5. [Performance](#performance)
6. [Architecture](#architecture)
7. [Preprocessing Pipeline](#preprocessing-pipeline)
8. [Soil Properties](#soil-properties)
9. [OSSL Data Sources](#ossl-data-sources)
10. [Installation](#installation)
11. [Quick Start](#quick-start)
12. [Shiny Interface](#shiny-interface)
13. [Versioning](#versioning)
14. [Data Citations](#data-citations)
15. [Software Citation](#software-citation)
16. [License](#license)

---

## Overview

**autoSpectra** is an R package and interactive Shiny application for predicting soil physical and chemical properties from diffuse reflectance spectra. It provides OSSL-trained model families for VisNIR and MIR: a faithful replication of the [Safanelli et al. (2025)](https://doi.org/10.1371/journal.pone.0296545) reference Cubist+PCA pipeline, a PCA-VAE family that keeps that pipeline and replaces only Cubist with an ensemble of supervised VAEs from the [soilVAE](https://CRAN.R-project.org/package=soilVAE) package, and the same PCA-VAE trained on single instrument windows for partial-range instruments:

| Model family | Algorithm | Range | Bands | Partial-range instruments? |
|---|---|---|---|---|
| **OSSL_VisNIR_Cubist** | Cubist + PCA(120) + SNV     | 400–2500 nm (reflectance) | 1 051 @ 2 nm | requires full coverage |
| **OSSL_MIR_Cubist**    | Cubist + PCA(120) + SNV     | 600–4000 cm⁻¹  | 1 701 @ 2 cm⁻¹ | requires full coverage |
| **OSSL_VisNIR_PCAVAE** | soilVAE ensemble on PCA(120) + SNV | 400–2500 nm (reflectance) | 1 051 @ 2 nm | requires full coverage |
| **OSSL_MIR_PCAVAE**    | soilVAE ensemble on PCA(120) + SNV | 600–4000 cm⁻¹ | 1 701 @ 2 cm⁻¹ | requires full coverage |
| **OSSL_VisNIR_PCAVAE_800** / **_1000** / **_1350** | as OSSL_VisNIR_PCAVAE, one instrument window | 800, 1000 or 1350–2500 nm | 851 / 751 / 576 @ 2 nm | yes, via `predict_soil_agnostic()` |
| **OSSL_MIR_PCAVAE_700_3800** | as OSSL_MIR_PCAVAE, one instrument window | 700–3800 cm⁻¹ | 1 551 @ 2 cm⁻¹ | yes, via `predict_soil_agnostic()` |
| **OSSL_NeoSpectra_Cubist** / **_PCAVAE** | as above, NeoSpectra handheld library | 1350–2550 nm | 601 @ 2 nm | NeoSpectra only |
| **KSSL_VisNIR_Cubist** / **_PCAVAE**, **KSSL_MIR_Cubist** / **_PCAVAE** | as the OSSL family, trained on the KSSL rows only | as OSSL | as OSSL | as OSSL |

All OSSL families are trained on the [Open Soil Spectral Library v1.2](https://docs.soilspectroscopy.org), which each user downloads from the official OSSL storage; autoSpectra never redistributes it. All families are scored on the same held-out samples, and the Cubist family is also validated with Safanelli's 10-fold protocol and with leave-one-dataset-out (see [Performance](#performance)).

> **Which family?** For a full-range instrument, use **PCA-VAE**: under Safanelli's 10-fold cross-validation, with the same folds, spectra and principal components as Cubist, it is significantly more accurate in 116 of the 141 models of the paper's S2 Table, ties in 22 and is worse in 3 (see [Performance](#performance)), and its prediction intervals are as well calibrated. **Cubist** remains the reference replication of the published OSSL methodology and needs no Python. When an OSSL-like instrument covers only part of the range (e.g. a Bruker MPA from 800 nm, a PerkinElmer from 1000 nm), use **`predict_soil_agnostic()`**: it sends each spectrum to the PCA-VAE trained on the widest window the spectrum covers, and keeps nearly the full-range accuracy. No OSSL family transfers to an instrument with a different optical response, such as the NeoSpectra handheld: for such an instrument, train a model on its own labelled spectra with `train_sensor_model()` (see [Your own instrument](#your-own-instrument-train_sensor_model)). Always check `predict_applicability()`.

---

## Models

### OSSL_VisNIR_Cubist / OSSL_MIR_Cubist — Safanelli 2025 replication

```
Algorithm  : Cubist (committees in {1, 5, 10, 15, 20})
Predictors : First 120 principal components (Safanelli, p. 10)
Preprocess : SNV (Safanelli ring-trial selection, p. 9); VisNIR kept as
             reflectance, 400-2500 nm, like the OSSL models
Tuning     : 5-fold inner CV on a 2 000-sample subset (Safanelli, p. 9)
Targets    : log1p transform for skewed properties (Safanelli, p. 10)
AD         : Q-statistic on PCA residuals and Hotelling T2 on the PC scores;
             each limit = 99th percentile of the training spectra (see
             Methodology for why not Jackson-Mudholkar / the F limit)
Intervals  : normalised conformal prediction (Safanelli eqs. 9-10), calibrated
             on held-out samples; any level from 50% to 99%, default 68%
```

### OSSL_VisNIR_PCAVAE / OSSL_MIR_PCAVAE — PCA + soilVAE ensemble

```
Predictors : the Cubist family's: SNV spectra -> first 120 principal components,
             fitted on the same training rows (same 80/20 split, seed 42)
Algorithm  : 5 supervised VAEs from the soilVAE package (>= 0.2.0), seeds 1-5,
             encoder 256-128, latent 32, decoder 128-256; prediction = their mean
Training   : each PC standardised; early stopping on the prediction loss of a
             grouped tenth of the training rows (after a warm-up)
Targets    : log1p transform for skewed properties (as Cubist)
AD         : Q-statistic and Hotelling T2, empirical 99th percentiles (as Cubist)
Intervals  : normalised conformal; the width follows the spread of the 5
             members; any level from 50% to 99%, default 68%
Needs      : soilVAE (>= 0.2.0) with its Python backend (soilVAE::vae_install())
```

Only the regression model differs from the Cubist family, so the comparison between the two is like for like.

### Instrument windows — PCA-VAE for partial-coverage spectra

```
Families   : OSSL_VisNIR_PCAVAE_800, _1000, _1350   (800/1000/1350-2500 nm)
             OSSL_MIR_PCAVAE_700_3800                 (700-3800 cm-1)
Pipeline   : the OSSL PCA-VAE above, restricted to the window: every OSSL
             spectrum that covers the window, same split seed, same PCA setup
Use        : predict_soil_agnostic(df, "visnir") / ("mir"): each spectrum goes
             to the widest window it fully covers (the full-range PCA-VAE
             first); the output names the model used in `model_family`
AD, PIs    : as the PCA-VAE family
```

The windows are those of the OSSL instruments (Foss from 400 nm, Bruker MPA from 800 nm, PerkinElmer from 1000 nm, NeoSpectra-like ranges from 1350 nm; MIR instruments from 700 cm⁻¹ or up to 3800 cm⁻¹). They replace the masked soilVAE of earlier versions, a single network for every range, which was clearly less accurate (mean R² on full spectra 0.58 VisNIR, 0.67 MIR, against 0.735 and 0.802 for the PCA-VAE).

---

## Methodology — what we copied from Safanelli (2025) and what we added

The package implements **every key step** of Safanelli et al. (2025) so the Cubist family is a direct replication. The PCA-VAE families then change only the regression model, keeping the input, validation, applicability domain and uncertainty identical.

| Step | Safanelli 2025 (Cubist) | autoSpectra Cubist family | autoSpectra PCA-VAE families |
|---|---|---|---|
| Spectral range | 400–2500 / 600–4000 (OSSL models) | same | same, plus single instrument windows |
| Preprocessing | SNV alone | SNV (VisNIR on reflectance¹) | same as Cubist |
| Predictor | First 120 PCs | First 120 PCs | same 120 PCs, standardised |
| Algorithm | Cubist (committees ∈ {1, 5, 10, 15, 20}) | same | ensemble of 5 supervised VAEs (soilVAE) |
| Inner tuning | 5-fold CV, ~2 000 sample subset | same | early stopping on a grouped tenth of the training rows |
| Validation | 10-fold CV with refitting; external ring-trial test | 10-fold CV with refitting; leave-one-dataset-out; held-out test | same as Cubist (same folds and split) |
| Target transform | log1p (skewed properties) | same | same |
| Applicability domain | Q-statistic, Jackson–Mudholkar UCL (99%) | Q-statistic and Hotelling T2, empirical 99th percentiles² | same as Cubist |
| Uncertainty | Normalised conformal, 68% | same, calibrated on held-out samples, 50–99% | same, width follows the spread of the 5 members |
| Repeat visits of a site (LUCAS 2009/2015) | not grouped | kept on one side of every split | same |
| **Partial-coverage instruments** | **no** (separate NeoSpectra model) | no | **yes** (one model per instrument window) |

¹ The paper does not state the VisNIR input; the OSSL model code uses reflectance from 400 nm. A 10-fold CV test on seven properties matched or beat absorbance 350–2500 nm on all but organic carbon (−0.003 cor²). The PCA can also be fitted, as in the OSSL code, on every library spectrum outside the test fold, with standardised bands (`train_ossl.R --pca-fit library --pca-scale`).
² On OSSL the Jackson–Mudholkar approximation flagged 11–40% of the *training* spectra as out of domain at a nominal 1%, so the limit is taken from the training distribution itself. The J–M value is still stored in each model's metrics file. The same holds for Hotelling T2: the F-distribution limit flagged 15–18% of the OSSL VisNIR library, because PC scores of soil spectra have heavy tails. With the empirical limit, a median of 1.0% of held-out spectra exceed it (299 models). The published models record this limit from v0.5.1 on (Zenodo [10.5281/zenodo.24513780](https://doi.org/10.5281/zenodo.24513780)); it was computed from each model's own training rows (`tools/t2_limits.R`). A model without it falls back to the Q-statistic alone.

The window families are what handles instruments that don't cover the full range (Bruker MPA 800–2500 nm, PerkinElmer 1000–2500 nm, etc.): `predict_soil_agnostic()` chooses among them for each spectrum.

---

## Supported Instruments

Spectra on any wavelength grid are resampled automatically to the model grid. That makes them *readable*, not necessarily *predictable*: the models work for instruments whose optical response resembles those in OSSL, and the window families also when such an instrument covers only part of the range. An instrument with a different response is a different matter; see [Transfer to a different instrument](#transfer-to-a-different-instrument-neospectra).

| Family | Partial coverage (e.g. NeoSpectra 1350-2550 nm) |
|---|---|
| `predict_soil_agnostic()` (PCA-VAE by window) | works — each spectrum goes to the widest window it covers (VisNIR from 400, 800, 1000 or 1350 nm; MIR 600–4000 or 700–3800 cm⁻¹) |
| `OSSL_VisNIR_Cubist` (Safanelli) | needs 350–2500 nm; missing bands are left out of SNV rather than invented, but accuracy drops |
| `OSSL_MIR_Cubist` (Safanelli) | needs 600–4000 cm⁻¹ |
| `OSSL_VisNIR_PCAVAE` / `OSSL_MIR_PCAVAE` | need the full range, like Cubist |
| your own family from `train_sensor_model()` | the range of your instrument |

### VisNIR instruments in OSSL v1.2

| Instrument | Manufacturer | Range (nm) |
|-----------|-------------|-----------|
| FieldSpec 3 / 4 | Malvern Panalytical (ASD) | 350–2500 |
| XDS Rapid Content Analyzer | Foss | 400–2500 |
| MPA FT-NIR | Bruker | 800–2500 |
| Spectrum One FT-NIR | PerkinElmer | 1000–2500 |
| MicroNIR | VIAVI Solutions | 908–1676 |
| Other diffuse-reflectance VisNIR | Various | ≥400 |

### MIR instruments in OSSL v1.2

| Instrument | Manufacturer | Technique | Range (cm⁻¹) |
|-----------|-------------|-----------|-------------|
| ALPHA FTIR | Bruker | ATR (diamond) | 400–4000 |
| Tensor 27 FTIR | Bruker | DRIFTS | 400–4000 |
| Spectrum Two FT-MIR | PerkinElmer | ATR | 400–4000 |
| Nicolet iS10 / iS50 | Thermo Fisher | ATR / DRIFTS | 400–4000 |
| Other FTIR instruments | Various | ATR / DRIFTS | ≥600 |

> **Bringing a new instrument?** `predict_soil()` resamples your spectra to the model grid. Then check `predict_applicability()`: it flags spectra unlike anything in the training data (Q-statistic on the PCA residuals and Hotelling T2 on the PC scores, for Cubist and PCA-VAE alike).

---

## Performance

Three evaluations, each answering a different question. The numbers below are generated from the trained models and the cross-validation runs, never typed in: `Rscript tools/s2_comparison.R` for the 10-fold comparison with Safanelli's S2 Table, `Rscript tools/benchmark_table.R` for the held-out, window, interval and leave-one-dataset-out tables.

| Question | Evaluation |
|---|---|
| Does the Cubist family faithfully replicate Safanelli et al. (2025), and is the PCA-VAE more accurate? | Random 10-fold CV with refitting, Safanelli's protocol, on all 141 models of the paper's S2 Table |
| How accurate are the instrument-window models, and are the intervals calibrated? | Held-out test set (same samples for every family) |
| How well do the models predict soils from a source they never saw? | Leave-one-dataset-out |

### 10-fold cross-validation: Safanelli et al. (2025), our Cubist and the PCA-VAE

Safanelli's protocol on both families: random 10-fold CV with refitting on OSSL v1.2 (the PCA, the Cubist committees and the networks are refitted without the test fold). Cubist and PCA-VAE use **the same folds**, spectra and principal components, so their ten fold R² values are paired; *p* is a two-sided exact Wilcoxon signed-rank test of them (0.002 is the smallest value ten folds allow).
`oc` and `n.tot` are on the log1p scale, as the paper reports them; pH and clay in original units.

| Sensor | Property | Safanelli S2 R² / RMSE / CCC | Our Cubist R² / RMSE / CCC | PCA-VAE R² / RMSE / CCC | ΔR² | *p* |
|---|---|---|---|---|---|---|
| VisNIR | `oc` | 0.94 / 0.24 / 0.97 | 0.944 / 0.24 / 0.971 | **0.965** / 0.19 / 0.982 | +0.020 | 0.002 |
| VisNIR | `n.tot` | 0.92 / 0.08 / 0.96 | 0.919 / 0.08 / 0.957 | **0.947** / 0.07 / 0.973 | +0.028 | 0.002 |
| VisNIR | `ph.h2o` | 0.90 / 0.43 / 0.95 | 0.898 / 0.43 / 0.946 | **0.936** / 0.34 / 0.967 | +0.037 | 0.002 |
| VisNIR | `clay.tot` | 0.82 / 6.61 / 0.89 | 0.813 / 6.65 / 0.893 | **0.872** / 5.51 / 0.930 | +0.059 | 0.002 |
| MIR | `oc` | 0.99 / 0.12 / 0.99 | 0.989 / 0.12 / 0.994 | **0.990** / 0.11 / 0.995 | +0.001 | 0.004 |
| MIR | `n.tot` | 0.97 / 0.06 / 0.98 | 0.969 / 0.06 / 0.984 | **0.971** / 0.05 / 0.985 | +0.002 | 0.002 |
| MIR | `ph.h2o` | 0.91 / 0.40 / 0.95 | 0.911 / 0.40 / 0.953 | **0.951** / 0.29 / 0.975 | +0.040 | 0.002 |
| MIR | `clay.tot` | 0.90 / 5.64 / 0.94 | 0.905 / 5.42 / 0.949 | **0.947** / 4.04 / 0.973 | +0.042 | 0.002 |

Over all 141 models of the S2 Table (every property of the OSSL VisNIR, MIR and NeoSpectra libraries and of the KSSL subsets):

| Library | Models | PCA-VAE better | Tie | PCA-VAE worse | Mean R² Safanelli | Mean R² our Cubist | Mean R² PCA-VAE |
|---|---|---|---|---|---|---|---|
| VisNIR, OSSL | 23 | 22 | 1 | 0 | 0.722 | 0.712 | 0.771 |
| MIR, OSSL | 42 | 41 | 1 | 0 | 0.831 | 0.828 | 0.882 |
| NeoSpectra, OSSL | 29 | 10 | 18 | 1 | 0.598 | 0.596 | 0.616 |
| VisNIR, KSSL | 6 | 6 | 0 | 0 | 0.862 | 0.859 | 0.894 |
| MIR, KSSL | 41 | 37 | 2 | 2 | 0.848 | 0.845 | 0.891 |
| **All** | **141** | **116** | **22** | **3** | | | |

Better / worse: Wilcoxon *p* < 0.05; tie otherwise.

- **The replication reproduces the paper**: median |ΔR²| between our Cubist and the S2 Table is 0.004 over the 141 models, and 122 of them are within 0.02. The largest gaps are small NeoSpectra sets (e.g. `na.ext`, 0.30 against 0.42). Our folds are a little stricter than the paper's, because repeat visits of a LUCAS point (2009/2015) are kept in the same fold.
- **The PCA-VAE is more accurate**, with everything but the regression model held fixed: significantly better in 116 of 141 models, by at least 0.01 R² in 113 and 0.03 in 82. The gains are largest for texture, pH, exchangeable cations and micronutrients. Where Cubist is already near 0.99 (MIR carbon and nitrogen) the gain is 0.001–0.002: consistent over folds but of no practical importance.
- **Small libraries are the exception**: on NeoSpectra (about 2,000 soils) 18 of 29 models are ties. The three losses are NeoSpectra Mehlich-3 P (−0.085, a property neither model predicts well) and KSSL MIR Mehlich-3 Ca and Mg (−0.011 each).

### Held-out test set (original units)

The published models. Cubist and PCA-VAE: the same 80/20 split, the same SNV spectra and the same 120 principal components (OSSL-style PCA: fitted on every library spectrum outside the test set, bands standardised); only the regression model differs.
LUCAS revisits of the same point are kept on one side of every split.

| Sensor | Property | PCA-VAE R² / RMSE / CCC | Cubist R² / RMSE / CCC |
|---|---|---|---|
| VisNIR | `oc` | **0.962** / 2.42 / 0.980 | 0.942 / 2.99 / 0.970 |
| VisNIR | `n.tot` | **0.937** / 0.13 / 0.967 | 0.904 / 0.17 / 0.948 |
| VisNIR | `ph.h2o` | **0.935** / 0.34 / 0.967 | 0.899 / 0.43 / 0.947 |
| VisNIR | `clay.tot` | **0.875** / 5.47 / 0.932 | 0.814 / 6.67 / 0.893 |
| MIR | `oc` | **0.989** / 1.39 / 0.995 | **0.989** / 1.38 / 0.995 |
| MIR | `n.tot` | **0.885** / 0.21 / 0.939 | 0.883 / 0.22 / 0.937 |
| MIR | `ph.h2o` | **0.950** / 0.30 / 0.974 | 0.908 / 0.40 / 0.951 |
| MIR | `clay.tot` | **0.945** / 4.08 / 0.972 | 0.900 / 5.51 / 0.946 |

On this single split the picture matches the 10-fold CV above: R² +0.04 to +0.06 for pH and clay, +0.02 to +0.03 for VisNIR carbon and nitrogen, and practically equal MIR organic carbon and nitrogen. The PCA-VAE intervals are as well calibrated as Cubist's (below).

### PCA-VAE by instrument window (R²)

Each window family is the PCA-VAE pipeline trained on one instrument window, on every OSSL spectrum that covers it; `predict_soil_agnostic()` sends each spectrum to the widest window it covers. Every window family uses the same split as the full-range model, so all of them are scored on the same held-out samples.

| VisNIR | 400–2500 nm | 800–2500 nm | 1000–2500 nm | 1350–2500 nm |
|---|---|---|---|---|
| `oc` | 0.962 | 0.961 | 0.961 | 0.955 |
| `n.tot` | 0.937 | 0.933 | 0.932 | 0.929 |
| `ph.h2o` | 0.935 | 0.930 | 0.929 | 0.924 |
| `clay.tot` | 0.875 | 0.860 | 0.855 | 0.831 |

| MIR | 600–4000 cm⁻¹ | 700–3800 cm⁻¹ |
|---|---|---|
| `oc` | 0.989 | 0.990 |
| `n.tot` | 0.885 | 0.886 |
| `ph.h2o` | 0.950 | 0.952 |
| `clay.tot` | 0.945 | 0.948 |

The window models keep almost all of the full-range accuracy: over all 23 VisNIR properties the mean R² is 0.735 for 400–2500 nm, 0.736 for 800–2500, 0.725 for 1000–2500 and 0.713 for 1350–2500 nm; over the 42 MIR properties it is 0.802 for both MIR windows (Cubist on full spectra: 0.665 VisNIR, 0.703 MIR). Losses grow as the window narrows, mostly for texture. These are OSSL spectra; a real instrument with a different optical response is tested below.

Total sulfur looked like a second loss on this held-out split (VisNIR `s.tot`: R² 0.80 full range, 0.60 from 1350 nm), but that gap is one sample: a KSSL soil with 4.6% S that the window models over-predict and flag as low confidence. Without it the 1350–2500 nm model scores 0.80. A random 10-fold CV of `s.tot` (19,807 samples, same PCA and ensemble settings) shows no loss at all:

| VisNIR `s.tot`, 10-fold CV | R² (mean ± SD over folds) | Worst fold | CCC | R² on log1p scale |
|---|---|---|---|---|
| 400–2500 nm | 0.680 ± 0.257 | 0.00 | 0.840 | 0.829 |
| 800–2500 nm | 0.703 ± 0.262 | −0.02 | 0.850 | 0.825 |
| 1000–2500 nm | 0.763 ± 0.115 | 0.50 | 0.864 | 0.829 |
| 1350–2500 nm | 0.759 ± 0.077 | 0.64 | 0.855 | 0.818 |

Original-scale R² for `s.tot` depends on whether the few high-sulfur soils fall in the test fold. That is a limit of the property, not of the windows; the narrow windows are, if anything, steadier.

### Prediction intervals and applicability domain (held-out rows)

| Sensor | Property | Family | 68% interval covers | 90% covers | 95% covers | Flagged out of domain |
|---|---|---|---|---|---|---|
| VisNIR | `oc` | PCA-VAE | 68.5% | 90.0% | 95.2% | 0.9% |
| VisNIR | `n.tot` | PCA-VAE | 67.9% | 90.3% | 95.6% | 1.1% |
| VisNIR | `ph.h2o` | PCA-VAE | 69.6% | 90.5% | 95.2% | 0.8% |
| VisNIR | `clay.tot` | PCA-VAE | 66.6% | 89.9% | 95.1% | 0.9% |
| VisNIR | `oc` | Cubist | 68.1% | 90.5% | 95.7% | 0.9% |
| VisNIR | `n.tot` | Cubist | 68.6% | 90.3% | 95.2% | 1.1% |
| VisNIR | `ph.h2o` | Cubist | 69.0% | 90.8% | 96.0% | 0.8% |
| VisNIR | `clay.tot` | Cubist | 67.8% | 89.7% | 94.4% | 0.9% |
| MIR | `oc` | PCA-VAE | 68.1% | 90.2% | 95.0% | 1.0% |
| MIR | `n.tot` | PCA-VAE | 68.9% | 90.3% | 95.0% | 1.0% |
| MIR | `ph.h2o` | PCA-VAE | 67.2% | 90.2% | 95.0% | 1.2% |
| MIR | `clay.tot` | PCA-VAE | 69.1% | 91.2% | 95.5% | 0.9% |
| MIR | `oc` | Cubist | 68.1% | 89.8% | 94.8% | 1.0% |
| MIR | `n.tot` | Cubist | 69.0% | 90.4% | 94.9% | 1.0% |
| MIR | `ph.h2o` | Cubist | 67.8% | 89.8% | 95.2% | 1.2% |
| MIR | `clay.tot` | Cubist | 68.0% | 90.4% | 95.2% | 0.9% |

### Cubist under leave-one-dataset-out (original units)

Each fold withholds one whole OSSL contributing dataset. R² and CCC are weighted by fold size, RMSE pooled. These Cubist figures come from an earlier run whose PCA was fitted on the training rows only; leave-one-dataset-out for every family (Cubist, PCA-VAE and the window models) and every property measured in at least two datasets, with the published settings, is running and will replace this table.

| Sensor | Property | R² | RMSE | CCC | R² median fold | Folds |
|---|---|---|---|---|---|---|
| VisNIR | `oc` | 0.631 | 8.86 | 0.752 | 0.597 | 4 |
| VisNIR | `n.tot` | 0.469 | 0.47 | 0.636 | 0.715 | 3 |
| VisNIR | `ph.h2o` | 0.631 | 0.83 | 0.807 | 0.325 | 4 |
| VisNIR | `clay.tot` | 0.412 | 12.43 | 0.634 | 0.387 | 4 |
| MIR | `oc` | 0.829 | 5.44 | 0.886 | 0.830 | 7 |
| MIR | `n.tot` | 0.682 | 0.36 | 0.797 | 0.500 | 8 |
| MIR | `ph.h2o` | 0.758 | 0.64 | 0.867 | 0.624 | 8 |
| MIR | `clay.tot` | 0.156 | 15.45 | 0.602 | 0.229 | 8 |

This is the hardest test, and the numbers fall accordingly: a model trained on OSSL is markedly less accurate on a contributing dataset it has never seen. Expect these figures, not the random-split ones, when your samples come from a region or laboratory absent from OSSL. Folds disagree strongly (see the median), so performance depends on how close your soils are to the training data; check `predict_applicability()`.

### Transfer to a different instrument (NeoSpectra)

The instrument-range results above crop OSSL spectra, so they test coverage, not a different instrument. To test a genuinely different sensor, the VisNIR models were applied, without any retraining, to the OSSL NeoSpectra library: 8,067 scans of 2,106 soils on 13 NeoSpectra handheld units (1350–2550 nm), none of them in the VisNIR training data. `predict_soil_agnostic()` sends every scan to `OSSL_VisNIR_PCAVAE_1350`, the only window they cover.

| `OSSL_VisNIR_PCAVAE_1350` on NeoSpectra scans | oc | n.tot | ph.h2o | clay.tot |
|---|---|---|---|---|
| R², all scans | −111.5 | −80.0 | −1.89 | −0.54 |
| Squared correlation, all scans | 0.19 | 0.07 | 0.11 | 0.20 |
| Flagged out of domain, Q-statistic | 62.9% | 62.1% | 95.1% | 89.6% |
| Flagged out of domain, Q or Hotelling T2 | 84.6% | 96.4% | 95.5% | 90.1% |
| R², scans inside the domain (Q only) | −109.8 | −49.8 | −0.92 | −0.84 |
| R², scans inside the domain (Q and T2) | −78.0 | −65.0 | −0.77 | −0.85 |

**The models do not transfer to this instrument, and the applicability domain does not catch every failure.** After the same preprocessing, NeoSpectra spectra have a different baseline and shallower absorption bands than the benchtop instruments in OSSL (correlation of the mean spectra 0.63). The Q-statistic alone flags 62–95% of the scans; adding Hotelling T2 raises this to 85–96% (organic carbon 85%, nitrogen 96%, pH 96%, clay 90%). The scans that pass both checks are still predicted badly (R² below 0 for all four properties). A spectrum inside the domain is therefore not proof that the model applies: both checks measure how far a spectrum is from the library's spectra, and a constant instrument offset can leave a spectrum close to them yet shift every prediction. Safanelli et al. (2025) also needed a separate model calibrated on NeoSpectra spectra (S2 Table; external test in their Table 6); `train_sensor_model()` builds one, see the next section.

### Your own instrument: `train_sensor_model()`

An instrument whose spectra differ from OSSL's needs a model calibrated on its own labelled spectra. `train_sensor_model()` builds one on the OSSL grid within the instrument's range, with the OSSL preprocessing, and saves it as a family that `predict_soil()`, `predict_applicability()` and the intervals use like any other.

To choose the model, the OSSL NeoSpectra library was used as a test bench: organic carbon, 1350–2500 nm, the 90 African soils of Safanelli's Table 6 and a random 20% of the US soils held out, every scan of a soil kept on one side, R² of the per-soil mean prediction in original units (USA / Africa). The rest of the US soils, or a random subset of them, were used for training.

| Training soils | Cubist + PCA | PCA-VAE, instrument only | PCA-VAE, OSSL pretraining + fine-tuning |
|---|---|---|---|
| 150 | **0.66 / 0.71** | 0.60 / 0.25 | 0.62 / 0.53 |
| 400 | **0.73 / 0.81** | 0.71 / −1.10 | 0.57 / 0.60 |
| 1613 | 0.78 / **0.78** | **0.82** / 0.77 | 0.78 / 0.78 |

- With many soils the PCA-VAE is the more accurate model; with few it is less accurate than Cubist and can fail on soils from another region (Africa, 400 soils). `model = "auto"` (default) therefore uses Cubist below 1000 distinct soils and the PCA-VAE from 1000 on.
- Pretraining on OSSL cut to the instrument's range and then fine-tuning on the instrument (`pretrain = TRUE`) softened the network's failures with few soils but did not beat Cubist there, nor the instrument-only network with many soils. It is available but off by default.
- None of the US or African NeoSpectra soils is in the OSSL VisNIR library (checked on their lab values), so pretraining saw none of the test soils.

```r
# spectra: Soil_ID + one column per wavelength; y: lab values in original units
train_sensor_model(spectra, y = lab$oc, prop = "oc", family_id = "MyNeoSpectra",
                   groups = soil_of_scan,         # the soil each scan belongs to
                   out_dir = "models")
predict_soil(new_spectra, "MyNeoSpectra", model_dir = "models", with_pi = TRUE)
```

### Safanelli's independent evaluation (Tables 4 and 5)

The paper also scored its fixed OSSL models on external ring-trial soils: 60 soils scanned on 12 VisNIR instruments (Table 4) and 70 soils on 20 MIR instruments (Table 5), reporting the mean over instruments. Those spectra are not part of OSSL, so this test cannot be reproduced here, and it is not the same as leave-one-dataset-out. The paper's R² (Eq. 6) is also defined differently from ours (1 − SSE/SST) and cannot go negative. For reference, the Cubist + OSSL means:

| Property | VisNIR R² / RMSE / CCC (12 instruments) | MIR R² / RMSE / CCC (20 instruments) |
|---|---|---|
| Organic carbon (%) | 0.78 / 5.27 / 0.82 | 0.96 / 1.16 / 0.95 |
| pH (H₂O) | 0.50 / 1.22 / 0.65 | 0.75 / 0.72 / 0.84 |
| Clay (%) | 0.25 / 13.35 / 0.38 | 0.62 / 8.82 / 0.74 |
| K exchangeable (cmolc/kg) | 0.23 / 1.09 / 0.28 | 0.49 / 0.44 / 0.61 |

### What changed in November 2026

An audit found and fixed several defects that the earlier metrics did not show, because they were scored on full-range spectra. Among them: a fake absorbance spike in VisNIR spectra starting above 350 nm (which had also broken one leave-one-dataset-out fold: VisNIR `n.tot` went from R² −4.34 to +0.47), a masked soilVAE (since replaced by the window models) that never saw at training the partial spectra it receives at prediction, an applicability-domain limit rounded to zero, and over-wide prediction intervals. See NEWS.md.

---

## Architecture

### PCA-VAE — ensemble of supervised VAEs on principal components

```
Input: 120 standardised principal components of the SNV spectrum
        │
  Encoder Dense 256 → 128
        │
  z (32-D latent)
        │
   ┌────┴──────────────────────┐
   ▼ Decoder 128 → 256 → 120   ▼ Prediction head → 1
   (reconstructs the PCs)      (log1p / z-scored target)
```

Five members (seeds 1–5) from the [soilVAE](https://CRAN.R-project.org/package=soilVAE) package; the prediction is their mean, and their spread sets the width of the conformal interval. The window families use the same network on the PCs of their window.

### Applicability Domain

Every family reports, per spectrum, a score, the model's limit, and whether the spectrum is inside it (`predict_applicability()`).

```
Q = squared residual of the spectrum after projection on the 120 PCs
in domain if Q <= 99th percentile of Q over the training spectra
(Cubist, PCA-VAE and the window families alike)
```

Out-of-domain predictions are not blocked, but they should not be trusted, whatever the prediction interval says.

PCA-VAE families also report how much the five ensemble members disagree (`ensemble_sd`) and flag the prediction as `confident = FALSE` above the 95th percentile of that spread on the calibration data. A spectrum can be in domain and still flagged: on NeoSpectra electrical conductivity, one test sample of 748 made 67% of the squared error in original units, and the most disputed 5% of samples held most of the largest errors. Treat flagged predictions of skewed properties (ec, extractable Na, Mn) with care.

---

## Preprocessing Pipeline

Following Safanelli (2025), SNV is the canonical preprocessing step — it was selected through the SS4GG ring-trial experiment because it removes both instrument-specific scaling and particle-size scattering (the two dominant sources of inter-instrument variation) without sacrificing predictive accuracy.

```
                        Cubist family (Safanelli replication)
   VisNIR (reflectance)                        MIR (absorbance)
   --------------------                        ----------------
   Raw reflectance R                           Raw absorbance A
           |                                            |
           v                                            |
   A = -ln(R)                                           |
           |                                            |
           v                                            v
        SNV (per-row)                              SNV (per-row)
           |                                            |
           v                                            v
   PCA (first 120 PCs) - same per family             same
           |                                            |
           v                                            v
   Cubist (committees from inner CV)                 same
```

Bands an instrument does not cover stay missing through every step and are set to 0 only at the end, so they never distort the absorbance conversion or the SNV. The PCA-VAE and window families use the same preprocessing as Cubist:

```r
# Cubist (Safanelli replication)
preprocess = c("ABSORBANCE", "SNV")          # VisNIR
preprocess = c("SNV")                         # MIR
```

A spectrum that covers only part of the range is not padded or imputed: `predict_soil_agnostic()` sends it to the window family whose whole grid it covers, so every model sees only complete spectra.

---

## Soil Properties

### OSSL Level-1 Harmonized Properties (34 targets)

| Category | Variable | Unit | Description |
|----------|----------|------|-------------|
| **Organic** | `oc` | w.pct | Organic carbon |
| | `c.tot` | w.pct | Total carbon |
| | `n.tot` | w.pct | Total nitrogen |
| **Texture** | `clay.tot` | w.pct | Total clay (<0.002 mm) |
| | `silt.tot` | w.pct | Total silt (0.002–0.05 mm) |
| | `sand.tot` | w.pct | Total sand (0.05–2 mm) |
| **Reaction** | `ph.h2o` | index | pH in water |
| | `ph.cacl2` | index | pH in CaCl₂ |
| | `caco3` | w.pct | Calcium carbonate |
| | `acidity` | cmolc/kg | Exchangeable acidity |
| | `ec` | dS/m | Electrical conductivity |
| **Physical** | `bd` | g/cm³ | Bulk density |
| | `aggstb` | w.pct | Aggregate stability |
| | `awc.33.1500kPa` | w.frac | Available water content |
| | `wr.33kPa` | w.pct | Water retention at 33 kPa |
| | `wr.1500kPa` | w.pct | Water retention at 1500 kPa |
| **Exchange** | `cec` | cmolc/kg | Cation exchange capacity |
| | `ca.ext` | cmolc/kg | Exchangeable calcium |
| | `k.ext` | cmolc/kg | Exchangeable potassium |
| | `mg.ext` | cmolc/kg | Exchangeable magnesium |
| | `na.ext` | cmolc/kg | Exchangeable sodium |
| **Micronutrients** | `p.ext` | mg/kg | Extractable phosphorus |
| | `fe.ext` | mg/kg | Extractable iron |
| | `fe.dith` | w.pct | Crystalline iron (dithionite) |
| | `fe.ox` | w.pct | Amorphous iron (oxalate) |
| | `al.ext` | cmolc/kg | Exchangeable aluminum |
| | `al.ox` | w.pct | Amorphous aluminum |
| | `al.dith` | w.pct | Crystalline aluminum |
| | `mn.ext` | mg/kg | Extractable manganese |
| | `zn.ext` | mg/kg | Extractable zinc |
| | `cu.ext` | mg/kg | Extractable copper |
| | `b.ext` | mg/kg | Extractable boron |
| | `s.tot` | w.pct | Total sulfur |
| | `s.ext` | mg/kg | Extractable sulfur |

---

## OSSL Data Sources

All OSSL v1.2 data are publicly available under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Artifacts distributed from this repository that are **derived from OSSL v1.2** —
the trained models published on Zenodo ([10.5281/zenodo.23950129](https://doi.org/10.5281/zenodo.23950129)) —
inherit that attribution requirement: please cite the OSSL as described in
[Data Citations](#data-citations).

| Component | File | URL |
|-----------|------|-----|
| VisNIR spectra (L0) | `ossl_visnir_L0_v1.2.csv.gz` | `https://storage.googleapis.com/soilspec4gg-public/ossl_visnir_L0_v1.2.csv.gz` |
| MIR spectra (L0) | `ossl_mir_L0_v1.2.csv.gz` | `https://storage.googleapis.com/soilspec4gg-public/ossl_mir_L0_v1.2.csv.gz` |
| Soil lab data (L1) | `ossl_soillab_L1_v1.2.csv.gz` | `https://storage.googleapis.com/soilspec4gg-public/ossl_soillab_L1_v1.2.csv.gz` |

```r
ossl_download()          # VisNIR + MIR + soillab (~2 GB)
ossl_download("visnir")  # VisNIR only
ossl_download("mir")     # MIR only
```

| OSSL Level | Description |
|-----------|-------------|
| **L0** | Original contributed spectra (used for training) |
| **L1** | Harmonized soil lab values in common units (training targets) |

---

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install autoSpectra
devtools::install_github("HugoMachadoRodrigues/autoSpectra")
```

### Core dependencies

| Package | Role |
|---------|------|
| `shiny`, `shinyWidgets` | Interactive interface |
| `ggplot2` | Visualization |
| `keras3` + `tensorflow` | soilVAE training and inference (optional; Cubist families do not need it) |
| `Cubist`, `irlba` | Cubist + PCA families |
| `httr`, `data.table` | OSSL data download and loading |
| `readxl`, `writexl` | Excel I/O |
| `MASS` | Mahalanobis distance (pseudoinverse) |

**Setting up Keras/TensorFlow:**

```r
install.packages("keras3")
keras3::install_keras()   # installs TensorFlow in a managed virtualenv
```

---

## Quick Start

### 1. Install and download pre-trained models

```r
# install.packages("remotes")
remotes::install_github("HugoMachadoRodrigues/autoSpectra")

library(autoSpectra)

# One-time download (VisNIR/MIR PCA-VAE + Cubist). Only the requested
# families are fetched, checked against their SHA256 sums.
download_ossl_models()

# Or pick exactly what you need
download_ossl_models(c("OSSL_VisNIR_PCAVAE", "OSSL_VisNIR_PCAVAE_800",
                       "OSSL_VisNIR_PCAVAE_1000", "OSSL_VisNIR_PCAVAE_1350"))  # VisNIR, any range
download_ossl_models("all")                                    # every published family
```

> The PCA-VAE, NeoSpectra and KSSL families are published from the models
> release after v0.4.1; until then, train them locally with `Rscript
> train_ossl.R` (PCA-VAE needs soilVAE >= 0.2.0 and its Python backend,
> `soilVAE::vae_install()`).

> Models are extracted under `default_model_dir()`, the package's folder in
> `tools::R_user_dir("autoSpectra", "data")` (override with
> `options(autoSpectra.model_dir = ...)`). VAE weights are loaded into the
> session cache on first prediction; Cubist bundles are read on demand from
> `<family>/models/<prop>_cubist.rds`.

### 2. Predict from new spectra

```r
library(autoSpectra)
df <- readxl::read_excel("my_visnir_spectra.xlsx")

# (a) Any instrument range: each spectrum goes to the PCA-VAE of the widest
#     window it covers (400, 800, 1000 or 1350-2500 nm); `model_family` says which
out_any <- predict_soil_agnostic(df, "visnir",
                                 properties = c("oc", "clay.tot", "ph.h2o"),
                                 with_pi = TRUE)

# (b) Safanelli 2025 Cubist baseline - same call, different family
#     Add with_pi = TRUE to also return 68% prediction intervals
out_cubist <- predict_soil(df, "OSSL_VisNIR_Cubist",
                           properties = c("oc", "clay.tot", "ph.h2o"),
                           with_pi    = TRUE)

# Applicability domain (Q-statistic and Hotelling T2 on the PCA)
app <- predict_applicability(df, "OSSL_VisNIR_Cubist", "oc")
plot_applicability(app)

# Visualisation
plot_spectra(df, family = get_family("OSSL_VisNIR_Cubist"))
plot_mean_spectrum(df, family = get_family("OSSL_VisNIR_Cubist"))
```

### 3. Pre-load (VAE families) for batch prediction

```r
preload_ossl_models("OSSL_VisNIR_PCAVAE")   # caches all VAE weights once
for (f in my_files) {
  df  <- readxl::read_excel(f)
  out <- predict_soil(df, "OSSL_VisNIR_PCAVAE")  # reuses in-memory cache
  writexl::write_xlsx(out, sub("\\.xlsx$", "_pred.xlsx", f))
}
```

### 4. Launch the interactive app

```r
run_autoSpectra()
```

### Advanced: retrain from OSSL source data

```r
ossl_download()                    # ~2 GB, cached locally
```

```bash
# Train the OSSL families (~6-10 h on CPU)
Rscript train_ossl.R

# Subsets:
Rscript train_ossl.R --only-pcavae           # both full-range PCA-VAE families
Rscript train_ossl.R OSSL_VisNIR_PCAVAE_800 --pca-fit library --pca-scale  # one window
Rscript train_ossl.R --only-cubist           # Safanelli baseline only
Rscript train_ossl.R OSSL_MIR_Cubist         # one family
Rscript train_ossl.R --cv10                  # random 10-fold CV with refitting (Safanelli's protocol)
Rscript train_ossl.R --cv10 --lodo           # leave-one-dataset-out
Rscript train_ossl.R --quick                 # 4 properties, 20 epochs, 3 folds

# README tables are generated from the results, never typed in
Rscript tools/benchmark_table.R --check --final-dir models_final
```

After training, package the outputs for Zenodo upload with:

```bash
Rscript upload_models.R
```

---

## Shiny Interface

The app provides a streamlined four-tab workflow:

| Tab | Description |
|-----|-------------|
| **Preview** | Sample count, detected bands, spectral range overlap check |
| **Spectrum Viewer** | Per-sample spectrum with model grid overlay |
| **Mean Spectrum** | Mean ± SD ribbon across all uploaded samples |
| **Predictions** | Predicted soil properties (table), applicability domain plot, Excel download |

**Model selection** is the first step: choose a family from the sidebar. "VisNIR / MIR - any instrument range" runs `predict_soil_agnostic()`, which picks the window model for each spectrum and adds a "Model (window)" column to the results. The badge shows the active algorithm. VAE weights are loaded into memory on first prediction; Cubist bundles are read on demand. The "Prediction intervals" checkbox enables Safanelli-style 68% PIs for the Cubist and PCA-VAE families.

**Supported upload formats**: `.xlsx`, `.xls`, `.csv`

**Expected file format**:

```
Soil_ID  | 350  | 352  | 354  | ... | 2500
---------|------|------|------|-----|------
Sample_1 | 0.12 | 0.13 | 0.14 | ... | 0.08
Sample_2 | 0.18 | 0.19 | 0.20 | ... | 0.11
```

Column headers must be the wavelength (nm) or wavenumber (cm⁻¹) positions. The `Soil_ID` column is configurable.

---

## Versioning

### Development version

- **Hotelling T2** joins the Q-statistic in `predict_applicability()` for the Cubist and PCA-VAE families: a spectrum is in domain when it passes both (empirical 99th-percentile limits)
- **Published models v0.5.1** (Zenodo [10.5281/zenodo.24513780](https://doi.org/10.5281/zenodo.24513780)): the v0.5.0 models with the T2 limit added to each metrics file; nothing else changed. Earlier model versions were removed from GitHub and remain archived on Zenodo

### v0.5.0 — 2026-11-23 *(current CRAN submission)*

**Sensor-agnostic prediction by instrument window + PCA-VAE families**

- **PCA-VAE families** (`OSSL_VisNIR_PCAVAE`, `OSSL_MIR_PCAVAE`, `OSSL_NeoSpectra_PCAVAE`): the Cubist pipeline (SNV, 120 PCs, same split) with Cubist replaced by a 5-seed soilVAE ensemble; normalised conformal intervals and a low-confidence flag from ensemble disagreement
- **Window families** (`OSSL_VisNIR_PCAVAE_800`, `_1000`, `_1350`, `OSSL_MIR_PCAVAE_700_3800`) and **`predict_soil_agnostic()`**, which routes each spectrum to the widest window it fully covers; they replace the masked soilVAE, which is no longer published
- **KSSL families** (`KSSL_VisNIR_Cubist`, `KSSL_MIR_Cubist`, `KSSL_VisNIR_PCAVAE`, `KSSL_MIR_PCAVAE`) and `train_sensor_model()` for your own instrument
- **Closer to the OSSL models**: reflectance 400–2500 nm + SNV for VisNIR, library PCA with scaled bands (`--pca-fit library --pca-scale`)
- **Published models v0.5.0** (Zenodo [10.5281/zenodo.24451470](https://doi.org/10.5281/zenodo.24451470)); `download_ossl_models()` fetches only the requested families from the GitHub release and checks their SHA256 sums
- **Redistribution guard**: Cubist 1.0's `$train` copy of the training data is stripped and audited; the release audit checks every published file
- Shiny app: "VisNIR / MIR – any instrument range" selector

### v0.4.0 — 2026-11-04

**Safanelli (2025) faithful redesign + masked soilVAE**

- **Four official families**: `OSSL_VisNIR`, `OSSL_MIR` (masked soilVAE) and `OSSL_VisNIR_Cubist`, `OSSL_MIR_Cubist` (Safanelli replication)
- **Cubist + PCA(120) + SNV** baseline that exactly replicates the OSSL methodology of Safanelli et al. (2025) — same preprocessing, 5-fold inner CV committee tuning, log1p target transform, Q-statistic UCL applicability domain, normalised conformal prediction at 68% PI
- **Masked soilVAE**: input shape `2 * d_spec` carrying spectrum + binary coverage mask; trained with `spectral_dropout` augmentation simulating the ranges of ASD / Foss / Bruker MPA / PerkinElmer / NeoSpectra — one model for every coverage pattern (range only; a different optical response is not covered, see Performance)
- **Spectral ranges aligned with Safanelli**: VisNIR 350-2500 nm (1 076 bands), MIR 600-4000 cm⁻¹ (1 701 bands)
- **10-fold CV with refitting** (Safanelli's outer loop) — `Rscript train_ossl.R --cv10` writes per-property `cv10_summary.json`
- **Unified `predict_applicability()`**: dispatches Q-statistic for Cubist or Mahalanobis for VAE, returns the same schema (`score`, `threshold`, `in_domain`)
- **Shiny app**: 4-family selector with algorithm badge, prediction-interval toggle for Cubist families
- New modules: `R/applicability.R`, `R/conformal.R`, `R/cubist_baseline.R`, `R/train_cv.R`
- **Audit fixes (November 2026)**: one preprocessing path for training and prediction; the masked soilVAE trained on the same partial spectra it is given at prediction; applicability-domain limits and conformal intervals calibrated (see NEWS.md for the full list)

### v0.3.0 — 2026-03-13

**OSSL-only sensor-agnostic redesign**

- **Two official models**: `OSSL_VisNIR` and `OSSL_MIR` — full OSSL v1.2 training, all instruments
- **Scientific sensor-agnostic methodology**: absorbance/SG-smooth/SNV preprocessing + soilVAE information bottleneck + latent applicability domain
- **10-fold CV stratified by OSSL dataset**: honest cross-instrument performance estimates (mean ± SD)
- **In-memory model cache** (`R/cache.R`): `get_cached_model()`, `preload_ossl_models()` — zero-overhead repeated prediction
- **Simplified Shiny UX**: two-model selector, lazy preload, applicability domain panel
- `train_ossl.R` rewritten: stratified k-fold CV, per-fold metrics, `metrics_summary.json`
- All R CMD check WARNINGs resolved; pkgdown site live
- Zenodo DOI registered: `10.5281/zenodo.19004686`

### v0.2.0 — 2026-03-13

**Major: R package + OSSL integration**

- Converted from standalone Shiny app to a proper R package
- Added OSSL v1.2 data download and integration
- Expanded to 34 OSSL L1 soil properties
- Two-step SG preprocessing (explicit smooth + derivative passes)
- MIR support, applicability domain, visualization functions
- roxygen2 documentation (44 exported functions)

### v0.1.0 — 2025

**Initial release: Shiny application**

- Interactive Shiny app for soil spectral prediction
- 23 soil properties, local training data
- soilVAE asymmetric autoencoder architecture

---

## Data Citations

If you use autoSpectra with OSSL data or the pre-trained models, please cite both the OSSL paper and the OSSL dataset:

### Open Soil Spectral Library (OSSL) — paper and Cubist methodology

```bibtex
@article{Safanelli2025OSSL,
  title   = {Open Soil Spectral Library (OSSL): Building reproducible soil
             calibration models through open development and community engagement},
  author  = {Safanelli, Jos{\'e} Lucas and Hengl, Tomislav and Parente, Leandro and
             Minarik, Robert and Bloom, Dellena E. and Todd-Brown, Katherine and
             Gholizadeh, Asa and de Sousa Mendes, Wanderson and Sanderman, Jonathan},
  journal = {PLOS ONE},
  year    = {2025},
  volume  = {20},
  number  = {1},
  pages   = {e0296545},
  doi     = {10.1371/journal.pone.0296545},
  url     = {https://docs.soilspectroscopy.org}
}
```

### OSSL v1.2 — training data

```bibtex
@dataset{Safanelli2021OSSLdata,
  title     = {Open Soil Spectral Library (training data and calibration models)},
  author    = {Safanelli, Jos{\'e} Lucas and Hengl, Tomislav and Sanderman, Jonathan and
               Parente, Leandro},
  year      = {2021},
  version   = {v1.2},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.7599269}
}
```

### TensorFlow / Keras

```bibtex
@software{tensorflow2015,
  title  = {{TensorFlow}: Large-Scale Machine Learning on Heterogeneous Systems},
  author = {Abadi, Mart{\'i}n and others},
  year   = {2015},
  url    = {https://www.tensorflow.org}
}
```

---

## Author

**Hugo Rodrigues**  
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8070--8126-A6CE39.svg?logo=orcid)](https://orcid.org/0000-0002-8070-8126)
[![Twitter](https://img.shields.io/badge/Twitter-@Hugo__MRodrigues-1DA1F2.svg?logo=twitter)](https://twitter.com/Hugo_MRodrigues)

---

## Software Citation

```bibtex
@software{Rodrigues2026autoSpectra,
  title   = {autoSpectra: Soil Spectral Modelling, Visualization and Prediction},
  author  = {Rodrigues, Hugo},
  year    = {2026},
  version = {0.5.0},
  doi     = {10.5281/zenodo.19004686},
  url     = {https://github.com/HugoMachadoRodrigues/autoSpectra},
  license = {MIT},
  note    = {ORCID: 0000-0002-8070-8126}
}
```

> Formatted citation also available via `citation("autoSpectra")` in R.

---

## Contributing

Contributions are welcome. To report bugs or request features, open an issue on GitHub.

---

## License

MIT © 2025–2026 Hugo Rodrigues. See [LICENSE](LICENSE).

OSSL data accessed by this package is distributed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

<p align="center">
  <img src="man/figures/logo.png" alt="autoSpectra logo" width="200"/>
</p>

<p align="center">
  Made with ❤️ and 🌱 soil science · <a href="https://twitter.com/Hugo_MRodrigues">@Hugo_MRodrigues</a>
</p>
