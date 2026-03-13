# Package index

## OSSL Data

Download and manage Open Soil Spectral Library data

- [`ossl_cache_dir()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_cache_dir.md)
  : Default local cache directory for OSSL data
- [`ossl_download()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_download.md)
  : Download all three OSSL components (VisNIR, MIR, soillab)
- [`ossl_download_file()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_download_file.md)
  : Download one OSSL component file to the cache directory
- [`ossl_join()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_join.md)
  : Join OSSL spectra with soil lab data for model training
- [`ossl_l1_labels`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_l1_labels.md)
  : Fancy display labels for OSSL L1 properties
- [`ossl_l1_properties`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_l1_properties.md)
  : All OSSL Level-1 harmonized soil property variable names (34
  targets)
- [`ossl_load_raw()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_load_raw.md)
  : Load a cached OSSL CSV (gzipped) as a data.table
- [`ossl_mir_instruments`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_mir_instruments.md)
  : MIR instruments contributing to OSSL v1.2
- [`ossl_mir_matrix()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_mir_matrix.md)
  : Extract OSSL MIR spectra as a matrix aligned to the OSSL standard
  grid
- [`ossl_prepare()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_prepare.md)
  : Convenience: load and join all OSSL data for a given sensor type
- [`ossl_soillab()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_soillab.md)
  : Extract OSSL Level-1 soil lab data for a set of properties
- [`ossl_visnir_instruments`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_visnir_instruments.md)
  : VisNIR instruments contributing to OSSL v1.2
- [`ossl_visnir_matrix()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_visnir_matrix.md)
  : Extract OSSL VisNIR spectra as a matrix aligned to the OSSL standard
  grid

## Pre-trained Models

Download pre-trained soilVAE models from Zenodo (one-time setup)

- [`download_ossl_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/download_ossl_models.md)
  : Download pre-trained OSSL soilVAE models from Zenodo
- [`models_available()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/models_available.md)
  : Check whether pre-trained models are available locally

## Model Training

Train soilVAE models on OSSL or custom datasets

- [`build_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/build_soilVAE.md)
  : Build the soilVAE asymmetric autoencoder model
- [`train_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_soilVAE.md)
  : Train a soilVAE model for a single soil property
- [`train_ossl_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_ossl_models.md)
  : Train soilVAE models for all properties of one model family using
  OSSL data
- [`load_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/load_soilVAE.md)
  : Load a saved soilVAE model and its scaler
- [`save_scaler()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/save_scaler.md)
  : Save a z-score scaler to RDS and JSON
- [`save_metrics_json()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/save_metrics_json.md)
  : Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles,
  latent stats)
- [`metrics_from_y()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/metrics_from_y.md)
  : Compute prediction metrics (RMSE, R², RPIQ)
- [`split_idx()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/split_idx.md)
  : Split indices into train / calibration / test sets

## In-Memory Model Cache

Load models once per session and reuse them without disk I/O

- [`get_cached_model()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/get_cached_model.md)
  : Retrieve a soilVAE model from the in-memory cache
- [`preload_ossl_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/preload_ossl_models.md)
  : Pre-load all soilVAE models for a family into memory
- [`list_cached_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/list_cached_models.md)
  : List models currently held in the in-memory cache
- [`clear_model_cache()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/clear_model_cache.md)
  : Clear the in-memory model cache

## Prediction

Predict soil properties and assess applicability domain

- [`predict_soil()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_soil.md)
  : Predict soil properties from a spectral data frame using soilVAE
  models
- [`predict_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_applicability.md)
  : Compute applicability-domain scores for new samples
- [`format_predictions()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/format_predictions.md)
  : Round predictions and rename columns to display labels
- [`get_family()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/get_family.md)
  : Lookup a model family from the registry
- [`family_matches()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/family_matches.md)
  : Test whether a family matches a sensor type
- [`property_label()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/property_label.md)
  : Get a display label for a soil property key

## Preprocessing

Spectral preprocessing functions

- [`preprocess_visnir()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/preprocess_visnir.md)
  : Canonical two-step SG preprocessing for VisNIR (reflectance input)
- [`preprocess_mir()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/preprocess_mir.md)
  : Canonical two-step SG preprocessing for MIR (absorbance input)
- [`reflectance_to_absorbance()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/reflectance_to_absorbance.md)
  : Convert reflectance to absorbance
- [`apply_pipeline()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_pipeline.md)
  : Apply a spectral preprocessing pipeline
- [`apply_sg_matrix()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_sg_matrix.md)
  : Apply Savitzky-Golay filter row-wise to a spectral matrix
- [`apply_snv()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_snv.md)
  : Apply Standard Normal Variate (SNV) row-wise
- [`parse_sg()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/parse_sg.md)
  : Parse SG parameter string of the form "SG(m,p)", "SG(m,p,d)", etc.
- [`sg_safe_params()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/sg_safe_params.md)
  : Compute safe SG parameters given available number of columns
- [`resample_to_grid()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/resample_to_grid.md)
  : Resample a spectral matrix to a target wavelength grid via linear
  interpolation
- [`get_wavelengths()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/get_wavelengths.md)
  : Extract wavelength/wavenumber positions from a data frame's column
  names

## Visualization

Plotting functions

- [`plot_spectra()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_spectra.md)
  : Plot one or more spectra from a data frame
- [`plot_mean_spectrum()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_mean_spectrum.md)
  : Plot a mean spectrum with ±1 SD ribbon from a spectral matrix
- [`plot_predictions()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_predictions.md)
  : Plot predicted vs observed values for a soil property
- [`plot_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_applicability.md)
  : Plot latent-space applicability domain scores

## Registry

Model registry and metadata

- [`model_registry`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/model_registry.md)
  : Official autoSpectra model family registry
- [`ossl_l1_properties`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_l1_properties.md)
  : All OSSL Level-1 harmonized soil property variable names (34
  targets)
- [`ossl_l1_labels`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_l1_labels.md)
  : Fancy display labels for OSSL L1 properties
- [`ossl_visnir_instruments`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_visnir_instruments.md)
  : VisNIR instruments contributing to OSSL v1.2
- [`ossl_mir_instruments`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_mir_instruments.md)
  : MIR instruments contributing to OSSL v1.2

## Shiny App

Launch the interactive application

- [`run_autoSpectra()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/run_autoSpectra.md)
  : Launch the autoSpectra Shiny application

## Utilities

Internal helper functions

- [`feature_minmax()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/feature_minmax.md)
  : Compute per-band min/max feature range from a preprocessed matrix
- [`dir_create()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dir_create.md)
  : Create a directory if it does not already exist
- [`vneed()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/vneed.md)
  : Safe shiny::validate wrapper
