## Resubmission

This is a resubmission. In this version I have:

* Fixed the test failures on r-devel-windows: `audit_release_archive()` built
  a regular expression from the temporary directory path, which broke when the
  Windows path contained a backslash followed by a digit. It now takes the
  relative paths from `list.files()` directly.
* Reworded "Pre-trained" in the Description. The other words flagged as
  possibly misspelled (OSSL, Safanelli, Zenodo, et al., nm, reflectance,
  variational autoencoders) are correct: a data library, an author name, a
  data repository, and technical terms.

## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

Checked locally (macOS, R 4.6.2) with `R CMD check --as-cran`, and on GitHub
Actions with R release on Linux, macOS and Windows and with R devel on Linux.

## Notes for the reviewer

* The pre-trained models (several hundred MB) are not part of the package.
  `download_ossl_models()` fetches them on request from Zenodo
  (doi:10.5281/zenodo.24451470) into `tools::R_user_dir("autoSpectra", "data")`;
  nothing is written elsewhere unless the user sets
  `options(autoSpectra.model_dir = ...)`. Examples that need these downloads,
  or a 'shiny' session, are wrapped in `\dontrun{}` or `if (interactive())`.
* The suggested packages 'keras3' and 'soilVAE' need a Python installation with
  TensorFlow; tests that use them are skipped when it is not available.
