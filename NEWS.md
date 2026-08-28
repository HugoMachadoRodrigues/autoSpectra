# autoSpectra 0.4.0

Redesign around the methodology of Safanelli et al. (2025), plus a set of
safeguards that stop OSSL training data from being redistributed inside
published model archives.

## Data provenance and redistribution (read this before republishing to Zenodo)

Model archives built by earlier versions of `upload_models.R` embedded the OSSL
training database inside the Cubist bundles, by two independent routes:

* `pca$scores` — the 120 PCA coordinates of every training sample
  (105,992 x 120 for MIR). Since the same file also carried `pca$rotation` and
  `pca$center`, `.reconstruct_pca()` inverted them back to the input spectra at
  100% of the retained variance.
* `model$data` — `Cubist::cubist()` retains a verbatim CSV copy of the entire
  training matrix, predictors *and* response, in this slot. For
  `OSSL_MIR_Cubist/oc` that was 307,275,015 bytes: one line per OSSL sample,
  carrying both the PC scores and the laboratory value in plain text.

OSSL v1.2 is CC-BY 4.0, so redistribution with attribution is permitted and no
licence was breached. The problem is that a file published as a *model* silently
contained a copy of the source database, and nothing in the pipeline objected.
Across 16 bundles this was 1,005.9 MB of payload where 12.7 MB was intended.

Neither slot is read at prediction time, so removing them changes nothing:
`.project_pca()` uses only `rotation`/`center`/`scale`, and Cubist predictions
are bit-identical without `$data` (verified on all 8 bundles).

Three layers now prevent this:

* **Write side** — `.strip_pca_scores()` and `.strip_cubist_data()` are applied
  to the response model and the conformal error model before `saveRDS()`.
  `pca$n_train` is kept as a provenance count.
* **Publish side** — `upload_models.R` now ships by allowlist
  (`.PUBLISHABLE_PATTERNS`) instead of by denylist, and *refuses* to build an
  archive containing an unrecognised file. The previous denylist excluded only
  `*.log` and `*.tmp`, matched nothing in the repository, and would have shipped
  the leak without complaint.
* **Verify side** — `audit_release_archive()` opens a built zip and inspects
  every artefact, failing on any file outside the allowlist, any array with a
  per-sample axis, or any character blob containing delimited numeric records.
  `upload_models.R` runs it automatically and deletes any archive that fails.

`ossl_download_file()` remains the only way OSSL data enters the package: it is
fetched at run time from the official bucket documented at
<https://docs.soilspectroscopy.org/db-access.html> into a user-local cache, and
is never bundled or re-hosted.

## New: four model families

Two algorithms across two spectral domains, all sharing the same preprocessing,
target transforms, validation strategy and applicability-domain machinery:

* `OSSL_VisNIR`, `OSSL_MIR` — masked soilVAE. A two-channel
  `[SNV spectrum | coverage mask]` input with a masked-MSE reconstruction loss,
  trained with `spectral_dropout()` augmentation that simulates partial-coverage
  instruments, so prediction does not depend on a specific band range.
* `OSSL_VisNIR_Cubist`, `OSSL_MIR_Cubist` — a faithful replication of the
  Safanelli et al. (2025) reference pipeline: SNV, PCA(120), Cubist with
  committee tuning by inner 5-fold CV, Q-statistic applicability domain
  (Jackson–Mudholkar, 99%), and normalised conformal prediction at 68%.

## Improvements

* Metrics are reported in the original units. Properties trained under `log1p`
  are back-transformed with `expm1()` before scoring, so RMSE, R², RPIQ and CCC
  are directly comparable with Safanelli's Tables 4 and 5.
* `.fit_pca()` uses `irlba::prcomp_irlba()` for large matrices — a randomised
  truncated SVD that is far cheaper in memory and time than a full `prcomp()`
  at rank 120 on a 130k x 1076 matrix.
* Training is resumable. A property counts as finished only when both its model
  and its `_metrics.json` exist, so an interrupted run continues from the next
  unfinished property; `--force` retrains regardless.
* `keras3::clear_session()` and `gc()` run between properties. Without them,
  accumulated backend state aborted the process while building the second model
  of a run (`condition_variable wait failed` on macOS).
* Spectral ranges follow Safanelli: VisNIR 350–2500 nm (1,076 bands),
  MIR 600–4000 cm⁻¹ (1,701 bands).

## Bug fixes

* `keras3::save_model()` is called with `overwrite = TRUE`. Without it a re-run
  trained each property for its full epoch budget and then discarded the result
  with "File already exists" — a failure that only surfaced after the expensive
  part.
* An interrupted OSSL download no longer masquerades as a valid cache.
  Downloads land on a `.part` path and are renamed only after the HTTP status
  and the gzip stream both verify. Note that R's `gzfile()` cannot detect this:
  it reads a truncated member up to the cut and reports a clean EOF, so the
  CRC32/ISIZE trailer is checked with `gzip -t`.
* `split_idx()` used sampled values as negative indices instead of positions,
  which overlapped the train and test sets.
* `upload_models.R` evaluated `normalizePath(MODEL_DIR)` after `setwd(MODEL_DIR)`,
  resolving a nonexistent `models/models`; the archive path prefix was stripped
  correctly only by accident.
* Migrated to `keras3` and the `.keras` model format.

## Testing

129 tests. New coverage for the redistribution guards, the release auditor
(including planted-leak red-team cases and the false-positive cases that a
size-based check got wrong), OSSL cache integrity, the save/overwrite/load
round-trip for models carrying custom `Sampling` and `masked_mse` objects, and
the PCA and conformal paths.

# autoSpectra 0.3.0

* OSSL-only sensor-agnostic redesign, CI/CD, full documentation, Zenodo
  deposit (DOI [10.5281/zenodo.19004686](https://doi.org/10.5281/zenodo.19004686)).

# autoSpectra 0.2.0

* R package structure with OSSL v1.2 integration.
