EC outlier sensitivity diagnostic

Successful run configuration:

- property: `ec`
- family: `OSSL_VisNIR_Masked_Aligned`
- CV: grouped stratified K-fold by `Soil_ID`
- folds: `3`
- epochs: `20`
- seed: `42`
- spectral treatment: remove top `0.5%` global spectral outliers
- target treatment: winsorize to `[0.5%, 99.5%]`
- mode: `cv_only`

Observed results from the successful run log:

- baseline rerun (same 3-fold / 20-epoch setup, untreated):
  - `R2 = -220.249`
  - `RMSE = 3.612`
  - `CCC = 0.087`
- trim + winsor:
  - `R2 = -12.490`
  - `RMSE = 1.071`
  - `CCC = 0.144`

Interpretation:

- The combined outlier treatment materially improved `ec` under the same diagnostic setup.
- However, `ec` remained a poor target after treatment, so outliers are a contributing factor, not a complete explanation.

Useful companion files:

- `results/vae_outlier_audit_target_summary.csv`
- `results/vae_outlier_audit_overlap.csv`
- `results/vae_outlier_audit_by_dataset.csv`
- `models/OSSL_VisNIR_Masked_Aligned/metrics_summary.json`
