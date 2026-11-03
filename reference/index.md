# Package index

## OSSL Data

Download the Open Soil Spectral Library from its official storage bucket
into a user-local cache. autoSpectra never redistributes OSSL data.

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
- [`ossl_l1_log_transform`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/ossl_l1_log_transform.md)
  : Properties that Safanelli log1p-transforms before fitting
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

Download the published model families from Zenodo (one-time setup)

- [`download_ossl_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/download_ossl_models.md)
  : Download pre-trained OSSL soilVAE models from Zenodo
- [`models_available()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/models_available.md)
  : Check whether pre-trained models are available locally
- [`audit_release_archive()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/audit_release_archive.md)
  : Audit a built release archive for redistributed OSSL training data

## Prediction

Predict soil properties, with intervals and applicability domain

- [`predict_soil()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_soil.md)
  : Predict soil properties from a spectral data frame
- [`predict_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_applicability.md)
  : Compute applicability-domain scores for new spectra
- [`format_predictions()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/format_predictions.md)
  : Round predictions and rename columns to display labels
- [`get_family()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/get_family.md)
  : Lookup a model family from the registry
- [`family_matches()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/family_matches.md)
  : Test whether a family matches a sensor type
- [`property_label()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/property_label.md)
  : Get a display label for a soil property key

## Preprocessing

Spectral preprocessing, resampling and coverage masks

- [`preprocess_visnir()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/preprocess_visnir.md)
  : Canonical Safanelli preprocessing for VisNIR (reflectance input)
- [`preprocess_mir()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/preprocess_mir.md)
  : Canonical Safanelli preprocessing for MIR (absorbance input)
- [`reflectance_to_absorbance()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/reflectance_to_absorbance.md)
  : Convert reflectance to absorbance
- [`apply_pipeline()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_pipeline.md)
  : Apply a spectral preprocessing pipeline
- [`apply_sg_matrix()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_sg_matrix.md)
  : Apply Savitzky-Golay filter row-wise to a spectral matrix
- [`apply_snv()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_snv.md)
  : Apply Standard Normal Variate (SNV) row-wise
- [`apply_snv_masked()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_snv_masked.md)
  : Apply SNV using only the masked (covered) bands of each row
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
- [`mask_from_coverage()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/mask_from_coverage.md)
  : Build a coverage mask from a resampled spectral matrix
- [`has_spectral_coverage()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/has_spectral_coverage.md)
  : Which rows carry any spectral signal at all

## Target Transforms

log1p transform for skewed properties and its bounded inverse, which
keeps predictions inside the calibration range

- [`apply_target_transform()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/apply_target_transform.md)
  : Forward log1p transform for properties whose policy asks for it
- [`invert_target_transform()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/invert_target_transform.md)
  : Invert the target transform without letting it explode
- [`target_bounds()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/target_bounds.md)
  : Training-scale bounds used to keep a back-transform from
  extrapolating
- [`metrics_original_scale()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/metrics_original_scale.md)
  : Score predictions in original units, with the bounded back-transform

## soilVAE Family

Masked, sensor-agnostic variational autoencoder

- [`build_maskedVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/build_maskedVAE.md)
  : Build a masked soilVAE for sensor-agnostic prediction
- [`train_maskedVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_maskedVAE.md)
  : Train a masked soilVAE for one soil property
- [`build_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/build_soilVAE.md)
  : Build a standard soilVAE
- [`train_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_soilVAE.md)
  : Train a (standard) soilVAE model for one soil property
- [`load_soilVAE()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/load_soilVAE.md)
  : Load a saved soilVAE model and its scaler
- [`spectral_dropout()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/spectral_dropout.md)
  : Random spectral-dropout augmentation for masked-VAE training
- [`train_ossl_models()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_ossl_models.md)
  : Train soilVAE models for all properties of one model family using
  OSSL data

## Cubist + PCA Family

Replication of Safanelli et al. (2025)

- [`train_cubist_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_cubist_pca.md)
  : Train one Cubist + PCA(120) model for a single soil property
- [`load_cubist_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/load_cubist_pca.md)
  : Load a saved Cubist + PCA bundle for a property
- [`predict_cubist_pca()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_cubist_pca.md)
  : Predict with a Cubist + PCA bundle

## Applicability Domain

Q-statistic (Cubist) and latent-space Mahalanobis distance (soilVAE)

- [`q_statistic()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/q_statistic.md)
  : Compute the Q-statistic for each row of a PCA residual matrix
- [`q_statistic_ucl()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/q_statistic_ucl.md)
  : Upper Confidence Limit for the Q-statistic (Jackson-Mudholkar 1979)
- [`predict_q_statistic()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_q_statistic.md)
  : Compute Q-statistic and in-domain flag for new spectra against a
  fitted PCA
- [`mahalanobis_d2()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/mahalanobis_d2.md)
  : Squared Mahalanobis distance against a multivariate Gaussian (mu,
  Sigma)
- [`mahalanobis_threshold()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/mahalanobis_threshold.md)
  : Chi-squared threshold for the Mahalanobis test in latent space

## Prediction Intervals

Conformal calibration

- [`train_normalized_conformal()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/train_normalized_conformal.md)
  : Train a normalised-conformal calibration for a Cubist + PCA model
- [`predict_with_pi()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/predict_with_pi.md)
  : Predict with prediction intervals using a normalised-conformal
  calibration
- [`absolute_conformal()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/absolute_conformal.md)
  : Simple absolute-residual conformal calibration (used by the VAE
  family)

## Validation

Cross-validation, including leave-one-dataset-out

- [`cv_train_property()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/cv_train_property.md)
  : Perform 10-fold CV with refitting for one (family, property)
  combination
- [`random_kfold()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/random_kfold.md)
  : Build random 10-fold indices (Safanelli default)
- [`stratified_kfold()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/stratified_kfold.md)
  : Build Leave-One-Dataset-Out (LODO) k-fold indices
- [`split_idx()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/split_idx.md)
  : Split indices into train / calibration / test sets
- [`metrics_from_y()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/metrics_from_y.md)
  : Compute prediction metrics (RMSE, Bias, R², RPIQ, CCC)

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

## Visualization

Plotting functions

- [`plot_spectra()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_spectra.md)
  : Plot one or more spectra from a data frame
- [`plot_mean_spectrum()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_mean_spectrum.md)
  : Plot a mean spectrum with +/-1 SD ribbon
- [`plot_predictions()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_predictions.md)
  : Plot predicted vs observed values for a soil property
- [`plot_applicability()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/plot_applicability.md)
  : Plot applicability-domain scores

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

Lower-level helpers used by the trainers

- [`save_scaler()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/save_scaler.md)
  : Save a z-score scaler to RDS and JSON
- [`save_metrics_json()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/save_metrics_json.md)
  : Save per-model metrics JSON (RMSE, R2, RPIQ, conformal quantiles,
  latent stats)
- [`feature_minmax()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/feature_minmax.md)
  : Compute per-band min/max feature range from a preprocessed matrix
- [`dir_create()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/dir_create.md)
  : Create a directory if it does not already exist
- [`vneed()`](https://HugoMachadoRodrigues.github.io/autoSpectra/reference/vneed.md)
  : Safe shiny::validate wrapper
