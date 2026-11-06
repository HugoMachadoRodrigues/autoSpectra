# autoSpectra <img src="man/figures/logo.png" align="right" height="138"/>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/autoSpectra)](https://CRAN.R-project.org/package=autoSpectra)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/autoSpectra)](https://CRAN.R-project.org/package=autoSpectra)
[![Version](https://img.shields.io/badge/version-0.4.0-blue.svg)](https://github.com/HugoMachadoRodrigues/autoSpectra/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19004686.svg)](https://doi.org/10.5281/zenodo.19004686)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3.svg?logo=R)](https://www.r-project.org/)
[![OSSL](https://img.shields.io/badge/data-OSSL%20v1.2-2e7d32.svg)](https://docs.soilspectroscopy.org)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8070--8126-A6CE39.svg?logo=orcid)](https://orcid.org/0000-0002-8070-8126)
[![R-CMD-check](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/R-CMD-check.yml)
[![pkgdown](https://github.com/HugoMachadoRodrigues/autoSpectra/actions/workflows/pkgdown.yml/badge.svg)](https://HugoMachadoRodrigues.github.io/autoSpectra/)
<!-- badges: end -->

> **Soil spectral modeling, visualization, and prediction** — models trained on the Open Soil Spectral Library (OSSL v1.2) for VisNIR and MIR spectra: a Cubist + PCA replication of Safanelli et al. (2025), a PCA + soilVAE ensemble that takes exactly the same input and is more accurate, and a masked soilVAE that also accepts instruments covering only part of the spectral range. `train_sensor_model()` builds the same kind of model for your own instrument.

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

**autoSpectra** is an R package and interactive Shiny application for predicting soil physical and chemical properties from diffuse reflectance spectra. It provides **six OSSL-trained model families** — three algorithms × two spectral domains: a faithful replication of the [Safanelli et al. (2025)](https://doi.org/10.1371/journal.pone.0296545) reference Cubist+PCA pipeline, a PCA-VAE family that keeps that pipeline and replaces only Cubist with an ensemble of supervised VAEs from the [soilVAE](https://CRAN.R-project.org/package=soilVAE) package, and a masked soilVAE for partial-range instruments:

| Model family | Algorithm | Range | Bands | Partial-range instruments? |
|---|---|---|---|---|
| **OSSL_VisNIR**        | masked soilVAE              | 350–2500 nm    | 1 076 @ 2 nm | yes (coverage mask) |
| **OSSL_MIR**           | masked soilVAE              | 600–4000 cm⁻¹  | 1 701 @ 2 cm⁻¹ | yes (coverage mask) |
| **OSSL_VisNIR_Cubist** | Cubist + PCA(120) + SNV     | 350–2500 nm    | 1 076 @ 2 nm | requires full coverage |
| **OSSL_MIR_Cubist**    | Cubist + PCA(120) + SNV     | 600–4000 cm⁻¹  | 1 701 @ 2 cm⁻¹ | requires full coverage |
| **OSSL_VisNIR_PCAVAE** | soilVAE ensemble on PCA(120) + SNV | 350–2500 nm | 1 076 @ 2 nm | requires full coverage |
| **OSSL_MIR_PCAVAE**    | soilVAE ensemble on PCA(120) + SNV | 600–4000 cm⁻¹ | 1 701 @ 2 cm⁻¹ | requires full coverage |

All six families are trained on the [Open Soil Spectral Library v1.2](https://docs.soilspectroscopy.org), which each user downloads from the official OSSL storage; autoSpectra never redistributes it. All families are scored on the same held-out samples, and the Cubist family is also validated with Safanelli's 10-fold protocol and with leave-one-dataset-out (see [Performance](#performance)).

> **Which family?** For a full-range instrument, use **PCA-VAE**: on the same split, spectra and principal components as Cubist it is more accurate for 6 of the 8 benchmark cases and ties on the other two (MIR organic carbon and nitrogen), and it returns prediction intervals. **Cubist** remains the reference replication of the published OSSL methodology and needs no Python. Use the **masked soilVAE** when an OSSL-like instrument covers only part of the range (e.g. a Foss XDS from 400 nm, a Bruker MPA from 800 nm): it is trained on simulated partial spectra and keeps its accuracy across those ranges, where Cubist degrades. No OSSL family transfers to an instrument with a different optical response, such as the NeoSpectra handheld: for such an instrument, train a model on its own labelled spectra with `train_sensor_model()` (see [Your own instrument](#your-own-instrument-train_sensor_model)). Always check `predict_applicability()`.

---

## Models

### OSSL_VisNIR — masked soilVAE (partial coverage, 350–2500 nm)

```
Family ID  : OSSL_VisNIR
Algorithm  : masked Variational Autoencoder, latent dim = 32
Domain     : Visible + Near-Infrared reflectance
Grid       : 350 – 2500 nm at 2 nm (1 076 bands)
Input      : [spectrum | binary coverage mask]   shape = 2 * 1076
Pipeline   : -ln(R) -> Savitzky-Golay smooth (window 11, order 2) -> SNV over covered bands
Training   : each sample also enters as 2 simulated partial-coverage copies
             (instrument starting at 350, 400, 800, 1000 or 1350 nm), each
             preprocessed exactly as predict_soil() would preprocess it
AD         : Mahalanobis distance in the 32-D latent space; limit = 95th
             percentile of the training samples
Intervals  : none (point predictions)
```

### OSSL_MIR — masked soilVAE (partial coverage, 600–4000 cm⁻¹)

```
Family ID  : OSSL_MIR
Algorithm  : masked Variational Autoencoder, latent dim = 32
Domain     : Mid-Infrared absorbance (ATR, DRIFTS)
Grid       : 600 - 4000 cm-1 at 2 cm-1 (1 701 bands)
Input      : [spectrum | binary coverage mask]   shape = 2 * 1701
Pipeline   : Savitzky-Golay smooth (window 11, order 2) -> SNV over covered bands
Training   : 2 simulated copies per sample (from 600, 650 or 700 cm-1, or up to 3800 cm-1)
AD         : Mahalanobis distance in the 32-D latent space (95th percentile limit)
Intervals  : none (point predictions)
```

### OSSL_VisNIR_Cubist / OSSL_MIR_Cubist — Safanelli 2025 replication

```
Algorithm  : Cubist (committees in {1, 5, 10, 15, 20})
Predictors : First 120 principal components (Safanelli, p. 10)
Preprocess : SNV (Safanelli ring-trial selection, p. 9); VisNIR converted to
             absorbance first
Tuning     : 5-fold inner CV on a 2 000-sample subset (Safanelli, p. 9)
Targets    : log1p transform for skewed properties (Safanelli, p. 10)
AD         : Q-statistic on PCA residuals; limit = 99th percentile of the
             training Q (see Methodology for why not Jackson-Mudholkar)
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
AD         : Q-statistic on PCA residuals, empirical 99th percentile (as Cubist)
Intervals  : normalised conformal; the width follows the spread of the 5
             members; any level from 50% to 99%, default 68%
Needs      : soilVAE (>= 0.2.0) with its Python backend (soilVAE::vae_install())
```

Only the regression model differs from the Cubist family, so the comparison between the two is like for like.

---

## Methodology — what we copied from Safanelli (2025) and what we added

The package implements **every key step** of Safanelli et al. (2025) so the Cubist family is a direct replication. The masked soilVAE family then changes only the predictive engine, keeping preprocessing, validation, applicability domain, and uncertainty as comparable as possible.

| Step | Safanelli 2025 (Cubist) | autoSpectra Cubist family | autoSpectra masked soilVAE family |
|---|---|---|---|
| Spectral range | 350–2500 / 600–4000 | same | same |
| Preprocessing | SNV alone | SNV (VisNIR: absorbance first¹) | SG smooth(11,2) → SNV over covered bands |
| Predictor | First 120 PCs | First 120 PCs | 32-D VAE latent |
| Algorithm | Cubist (committees ∈ {1, 5, 10, 15, 20}) | same | masked VAE |
| Inner tuning | 5-fold CV, ~2 000 sample subset | same | early stopping on held-out samples |
| Validation | 10-fold CV with refitting; external ring-trial test | 10-fold CV with refitting; leave-one-dataset-out; held-out test | held-out test, also per instrument range |
| Target transform | log1p (skewed properties) | same | same |
| Applicability domain | Q-statistic, Jackson–Mudholkar UCL (99%) | Q-statistic, empirical 99th percentile² | Mahalanobis on latent, empirical 95th percentile |
| Uncertainty | Normalised conformal, 68% | same, calibrated on held-out samples, 50–99% | none |
| Repeat visits of a site (LUCAS 2009/2015) | not grouped | kept on one side of every split | same |
| **Partial-coverage instruments** | **no** (separate NeoSpectra model) | no | **yes** (mask + simulated coverage) |

¹ The paper does not state an absorbance step for VisNIR; OSSL VisNIR is stored as reflectance, and SNV is applied to −ln(R).
² On OSSL the Jackson–Mudholkar approximation flagged 11–40% of the *training* spectra as out of domain at a nominal 1%, so the limit is taken from the training distribution itself. The J–M value is still stored in each model's metrics file.

The masked soilVAE is the **only family that handles instruments that don't cover the full range** (NeoSpectra 1350–2550 nm, Bruker MPA 800–2500 nm, etc.) without requiring a separate model.

---

## Supported Instruments

Spectra on any wavelength grid are resampled automatically to the model grid. That makes them *readable*, not necessarily *predictable*: the models work for instruments whose optical response resembles those in OSSL, and the masked soilVAE also when such an instrument covers only part of the range. An instrument with a different response is a different matter; see [Transfer to a different instrument](#transfer-to-a-different-instrument-neospectra).

| Family | Partial coverage (e.g. NeoSpectra 1350-2550 nm) |
|---|---|
| `OSSL_VisNIR` (masked VAE) | works — coverage encoded in mask channel |
| `OSSL_MIR` (masked VAE) | works |
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

> **Bringing a new instrument?** `predict_soil()` resamples your spectra to the model grid. Then check `predict_applicability()`: it flags spectra unlike anything in the training data (Q-statistic for Cubist, latent Mahalanobis distance for soilVAE).

---

## Performance

Three evaluations, each answering a different question. All numbers below are generated from the trained models by `Rscript tools/benchmark_table.R`, never typed in.

| Question | Evaluation |
|---|---|
| Does the Cubist family faithfully replicate Safanelli et al. (2025)? | Random 10-fold CV, compared with the paper's S2 Table |
| How do the families compare, and how robust is the masked soilVAE to partial coverage? | Held-out test set, also cut to the range of real instruments |
| How well do the models predict soils from a source they never saw? | Leave-one-dataset-out |

### Cubist vs Safanelli et al. (2025): 10-fold cross-validation

Same protocol on both sides: Cubist + PCA(120) + SNV, random 10-fold CV on OSSL v1.2.
`oc` and `n.tot` are on the log1p scale, as the paper reports them; pH and clay in original units.

| Sensor | Property | Scale | autoSpectra R² / RMSE / CCC | Safanelli S2 R² / RMSE / CCC |
|---|---|---|---|---|
| VisNIR | `oc` | log1p | 0.942 / 0.24 / 0.970 | 0.94 / 0.24 / 0.97 |
| VisNIR | `n.tot` | log1p | 0.923 / 0.08 / 0.959 | 0.92 / 0.08 / 0.96 |
| VisNIR | `ph.h2o` | original | 0.887 / 0.45 / 0.939 | 0.90 / 0.43 / 0.95 |
| VisNIR | `clay.tot` | original | 0.815 / 6.62 / 0.894 | 0.82 / 6.61 / 0.89 |
| MIR | `oc` | log1p | 0.989 / 0.12 / 0.994 | 0.99 / 0.12 / 0.99 |
| MIR | `n.tot` | log1p | 0.969 / 0.06 / 0.984 | 0.97 / 0.06 / 0.98 |
| MIR | `ph.h2o` | original | 0.912 / 0.39 / 0.953 | 0.91 / 0.40 / 0.95 |
| MIR | `clay.tot` | original | 0.909 / 5.32 / 0.951 | 0.90 / 5.64 / 0.94 |

The replication reproduces the paper: R² within 0.01–0.02, RMSE and CCC within rounding. pH on VisNIR is slightly lower, and our folds are a little stricter than the paper's, because repeat visits of a LUCAS point (2009/2015) are kept in the same fold.

### Held-out test set, all families (original units)

Cubist and PCA-VAE: the same 80/20 split, the same SNV spectra and the same 120 principal components; only the regression model differs. Masked soilVAE: 60/20/20 train/calibration/test with the same seed, so its test rows largely overlap.
LUCAS revisits of the same point are kept on one side of every split.

| Sensor | Property | PCA-VAE R² / RMSE / CCC | Cubist R² / RMSE / CCC | masked soilVAE R² / RMSE / CCC |
|---|---|---|---|---|
| VisNIR | `oc` | **0.954** / 2.66 / 0.976 | 0.939 / 3.06 / 0.969 | 0.866 / 4.55 / 0.914 |
| VisNIR | `n.tot` | **0.927** / 0.14 / 0.962 | 0.912 / 0.16 / 0.953 | 0.889 / 0.18 / 0.938 |
| VisNIR | `ph.h2o` | **0.935** / 0.35 / 0.966 | 0.887 / 0.45 / 0.940 | 0.849 / 0.53 / 0.917 |
| VisNIR | `clay.tot` | **0.868** / 5.61 / 0.929 | 0.821 / 6.55 / 0.896 | 0.804 / 6.84 / 0.885 |
| MIR | `oc` | 0.989 / 1.41 / 0.995 | **0.990** / 1.31 / 0.995 | 0.814 / 5.79 / 0.866 |
| MIR | `n.tot` | **0.885** / 0.22 / 0.939 | 0.882 / 0.22 / 0.937 | 0.856 / 0.24 / 0.912 |
| MIR | `ph.h2o` | **0.949** / 0.30 / 0.974 | 0.907 / 0.40 / 0.950 | 0.924 / 0.37 / 0.959 |
| MIR | `clay.tot` | **0.945** / 4.08 / 0.972 | 0.905 / 5.38 / 0.949 | 0.899 / 5.53 / 0.944 |

With the same input, the soilVAE ensemble is more accurate than Cubist in 6 of the 8 cases: R² +0.04 to +0.05 for pH and clay, +0.015 for VisNIR carbon and nitrogen. MIR nitrogen (+0.003) and MIR organic carbon (−0.001) are ties. Its intervals are as well calibrated as Cubist's (below). The masked soilVAE, which also has to handle partial spectra, remains the least accurate on full spectra: it underestimates carbon-rich soils, most visibly on MIR (on the log1p scale its R² there is 0.967), and its results vary more between training runs.

### soilVAE by instrument range (R², same test rows)

Each test spectrum is cut to the range of a real instrument and predicted through `predict_soil()`'s preprocessing.

| VisNIR | full 350 2500 | from 400 | from 800 | from 1000 | from 1350 |
|---|---|---|---|---|---|
| `oc` | 0.87 | 0.84 | 0.85 | 0.85 | 0.87 |
| `n.tot` | 0.89 | 0.88 | 0.87 | 0.88 | 0.87 |
| `ph.h2o` | 0.85 | 0.85 | 0.85 | 0.84 | 0.83 |
| `clay.tot` | 0.80 | 0.80 | 0.78 | 0.77 | 0.77 |

| MIR | full 600 4000 | from 650 | from 700 | to 3800 |
|---|---|---|---|---|
| `oc` | 0.81 | 0.82 | 0.82 | 0.82 |
| `n.tot` | 0.86 | 0.86 | 0.86 | 0.86 |
| `ph.h2o` | 0.92 | 0.92 | 0.92 | 0.92 |
| `clay.tot` | 0.90 | 0.90 | 0.90 | 0.90 |

This is what the soilVAE is for: its accuracy holds when an OSSL instrument's spectrum is cut to a shorter range, down to 1350–2500 nm. Cubist has no such mechanism; use it with full-range spectra. These are cropped OSSL spectra; a real instrument with a different optical response is tested below.

### Prediction intervals and applicability domain (held-out rows)

| Sensor | Property | Family | 68% interval covers | 90% covers | 95% covers | Flagged out of domain |
|---|---|---|---|---|---|---|
| VisNIR | `oc` | PCA-VAE | 67.6% | 90.5% | 95.2% | 1.1% |
| VisNIR | `n.tot` | PCA-VAE | 67.7% | 89.8% | 94.9% | 0.9% |
| VisNIR | `ph.h2o` | PCA-VAE | 69.8% | 89.9% | 95.0% | 0.9% |
| VisNIR | `clay.tot` | PCA-VAE | 65.8% | 89.3% | 94.2% | 1.2% |
| VisNIR | `oc` | Cubist | 67.5% | 90.4% | 95.3% | 1.1% |
| VisNIR | `n.tot` | Cubist | 67.9% | 89.2% | 94.9% | 0.9% |
| VisNIR | `ph.h2o` | Cubist | 68.9% | 90.8% | 95.4% | 0.9% |
| VisNIR | `clay.tot` | Cubist | 67.2% | 91.0% | 94.8% | 1.2% |
| MIR | `oc` | PCA-VAE | 67.5% | 89.9% | 94.9% | 0.9% |
| MIR | `n.tot` | PCA-VAE | 68.6% | 90.1% | 95.0% | 1.0% |
| MIR | `ph.h2o` | PCA-VAE | 68.3% | 90.3% | 95.0% | 1.2% |
| MIR | `clay.tot` | PCA-VAE | 70.4% | 90.9% | 95.1% | 0.8% |
| MIR | `oc` | Cubist | 67.8% | 90.6% | 95.0% | 0.9% |
| MIR | `n.tot` | Cubist | 68.3% | 90.8% | 95.0% | 1.0% |
| MIR | `ph.h2o` | Cubist | 68.0% | 90.8% | 95.1% | 1.2% |
| MIR | `clay.tot` | Cubist | 68.6% | 90.5% | 95.2% | 0.8% |

### Cubist under leave-one-dataset-out (original units)

Each fold withholds one whole OSSL contributing dataset. R² and CCC are weighted by fold size, RMSE pooled.

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

The instrument-range results above crop OSSL spectra, so they test coverage, not a different instrument. To test a genuinely different sensor, the VisNIR models were applied, without any retraining, to the OSSL NeoSpectra library: 8,067 scans of 2,106 soils on 13 NeoSpectra handheld units (1350–2550 nm), none of them in the VisNIR training data.

| VisNIR model on NeoSpectra scans | oc | n.tot | ph.h2o | clay.tot |
|---|---|---|---|---|
| soilVAE R² | −1.16 | −0.79 | 0.30 | −2.62 |
| Flagged out of domain | 100% | | | |

**The models do not transfer to this instrument.** After the same preprocessing, NeoSpectra spectra have a different baseline and shallower absorption bands than the benchtop instruments in OSSL (correlation of the mean spectra 0.63), and no augmentation used in training covers that. The applicability domain did its job and flagged every scan, so `predict_applicability()` warns rather than silently returning wrong values. Safanelli et al. (2025) also needed a separate model calibrated on NeoSpectra spectra (S2 Table; external test in their Table 6); `train_sensor_model()` builds one, see the next section.

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

An audit found and fixed several defects that the earlier metrics did not show, because they were scored on full-range spectra. Among them: a fake absorbance spike in VisNIR spectra starting above 350 nm (which had also broken one leave-one-dataset-out fold: VisNIR `n.tot` went from R² −4.34 to +0.47), a masked soilVAE that never saw at training the partial spectra it receives at prediction, an applicability-domain limit rounded to zero, and over-wide prediction intervals. See NEWS.md.

---

## Architecture

### soilVAE — masked variational autoencoder

```
Input [spectrum | coverage mask]  (2 x d bands)
        │
  Dense 512 → 256 → 128 → 64   (ReLU, batch norm, dropout)
        │
  z_mean, z_log_var  (32-D) ── sampling during training; z_mean at prediction
        │
   ┌────┴─────────────────────────┐
   ▼ Reconstruction head          ▼ Prediction head
 Dense 64 → 128 → 256 → d       Dense 256 → 128 → 64 → 1
 loss: MSE on covered bands     loss: MSE on the (log1p / z-scored) target
 weight 0.1                     weight 1.0
```

| Hyperparameter | Value |
|---------------|-------|
| Latent dimension | 32 |
| Optimizer | Adam, learning rate 1×10⁻³ |
| Max epochs / batch size | 60 / 128 |
| KL weight | annealed from 0 to 5×10⁻⁴ over the first 20 epochs |
| Early stopping | monitors the prediction loss on held-out samples (15% of training samples, never copies of fitted ones); starts after the KL warm-up; patience 12 |
| LR reduction | ×0.5 when the prediction loss plateaus (patience 5), min 1×10⁻⁵ |

### Applicability Domain

Both families report, per spectrum, a score, the model's limit, and whether the spectrum is inside it (`predict_applicability()`).

```
Cubist : Q = squared residual of the spectrum after projection on the 120 PCs
         in domain if Q <= 99th percentile of Q over the training spectra
soilVAE: D² = (z − μ)ᵀ Σ⁻¹ (z − μ) in the 32-D latent space
         in domain if D² <= 95th percentile of D² over the training spectra
```

Out-of-domain predictions are not blocked, but they should not be trusted, whatever the prediction interval says.

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

Bands an instrument does not cover stay missing through every step and are set to 0 only at the end, so they never distort the absorbance conversion, the smoothing or the SNV. The masked soilVAE family additionally Savitzky-Golay smooths the spectrum before SNV, and computes the SNV mean and SD from the covered bands only:

```r
# Cubist (Safanelli replication)
preprocess = c("ABSORBANCE", "SNV")          # VisNIR
preprocess = c("SNV")                         # MIR

# Masked soilVAE (partial coverage)
preprocess = c("ABSORBANCE", "SG_SMOOTH(11,2)", "SNV_MASKED")  # VisNIR
preprocess = c("SG_SMOOTH(11,2)", "SNV_MASKED")                # MIR
```

For the masked-VAE family, training-time augmentation (`spectral_dropout`) adds two copies of each sample restricted to the range of a real instrument (Foss from 400 nm, Bruker MPA from 800 nm, PerkinElmer from 1000 nm, NeoSpectra from 1350 nm, ...). Each copy goes through the same preprocessing call as `predict_soil()`, so it is exactly the input the model would receive from that instrument.

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

# One-time download from Zenodo (all four families)
download_ossl_models()

# Or pick exactly what you need
download_ossl_models(c("OSSL_VisNIR", "OSSL_VisNIR_Cubist"))   # VisNIR pair
```

> The PCA-VAE families (`OSSL_VisNIR_PCAVAE`, `OSSL_MIR_PCAVAE`) are not in
> the Zenodo record yet; until they are, train them locally with
> `Rscript train_ossl.R --only-pcavae` (needs soilVAE >= 0.2.0 and its Python
> backend, `soilVAE::vae_install()`).

> Models are extracted under `models/` (override with
> `options(autoSpectra.model_dir = ...)`). VAE weights are loaded into the
> session cache on first prediction; Cubist bundles are read on demand from
> `<family>/models/<prop>_cubist.rds`.

### 2. Predict from new spectra

```r
library(autoSpectra)
df <- readxl::read_excel("my_visnir_spectra.xlsx")

# (a) Masked VAE - also accepts spectra covering only part of 350-2500 nm
out_vae <- predict_soil(df, "OSSL_VisNIR",
                        properties = c("oc", "clay.tot", "ph.h2o"))

# (b) Safanelli 2025 Cubist baseline - same call, different family
#     Add with_pi = TRUE to also return 68% prediction intervals
out_cubist <- predict_soil(df, "OSSL_VisNIR_Cubist",
                           properties = c("oc", "clay.tot", "ph.h2o"),
                           with_pi    = TRUE)

# Applicability domain (auto-dispatches Q-statistic vs Mahalanobis)
app <- predict_applicability(df, "OSSL_VisNIR_Cubist", "oc")
plot_applicability(app)

# Visualisation
plot_spectra(df, family = get_family("OSSL_VisNIR"))
plot_mean_spectrum(df, family = get_family("OSSL_VisNIR"))
```

### 3. Pre-load (VAE families) for batch prediction

```r
preload_ossl_models("OSSL_VisNIR")          # caches all VAE weights once
for (f in my_files) {
  df  <- readxl::read_excel(f)
  out <- predict_soil(df, "OSSL_VisNIR")    # reuses in-memory cache
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
# Train all four families (GPU recommended for VAE; ~6-10 h on CPU)
Rscript train_ossl.R

# Subsets:
Rscript train_ossl.R --only-vae              # both masked VAEs only
Rscript train_ossl.R --only-cubist           # Safanelli baseline only
Rscript train_ossl.R OSSL_MIR_Cubist         # one family
Rscript train_ossl.R --cv10                  # random 10-fold CV with refitting (Safanelli's protocol)
Rscript train_ossl.R --cv10 --lodo           # leave-one-dataset-out
Rscript train_ossl.R --quick                 # 4 properties, 20 epochs, 3 folds

# README tables are generated from the results, never typed in
Rscript tools/benchmark_table.R --check
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

**Model selection** is the first step: choose one of the four families (VisNIR / MIR × masked soilVAE / Cubist+PCA) from the sidebar. The badge shows the active algorithm. VAE weights are loaded into memory on first prediction; Cubist bundles are read on demand. The "Prediction intervals" checkbox enables Safanelli-style 68% PIs for the Cubist families.

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

### v0.4.0 — 2026-11-04 *(current)*

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
  version = {0.4.0},
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
