# Temporal evaluation (Fig. S2B)

This directory contains the adopted temporal evaluation: early 2023–2024 development and late 2025–2026 evaluation. It is separate from the main 57-bag OOF analysis.

| File | Contents |
|---|---|
| `input.rds` | 1,512 model-input spectra; saved early/late partition; five early OOF folds; 1,010 early OOF and 502 late spectrum probabilities |
| `preprocessing.rds` | Saved mean imputation, zero-variance removals, normalization and five-component PCA parameters |
| `model.json` | Adopted final XGBoost model in native JSON format |
| `configuration.json` | Model settings, seed and recorded R/package versions |
| `methods.R` | Preprocessing, model, top-3 aggregation and early-only threshold-selection functions |
| `predict_saved.R` | Reloaded-model inference and comparison with the supplied predictions and bag scores |

The 22 late bag scores and fixed threshold **0.98162150382995605** remain in `supplement/data/s2b_temporal_late_bags.csv` and `s2b_temporal_threshold.csv`. The plotting code remains `supplement/code/plot_s2b.py`.

## Reuse the saved model

Use R and xgboost versions recorded in `configuration.json`, plus jsonlite. From the repository root, choose a new output directory:

```sh
Rscript --vanilla supplement/analysis/s2b/predict_saved.R . checked_s2b
```

This uses the saved numerical preprocessing parameters and native model; it does not call `prep`, train a model or select a threshold. The packaged export reproduced all 502 saved late probabilities exactly and all 22 bag classifications. Bag scores agree within 1e-14 (CSV serialization precision). This verifies reuse of this model, not numerical identity with an earlier model run or stability under other software settings.

## Read the inputs

```r
d <- readRDS("supplement/analysis/s2b/input.rds")
head(d$metadata)
dim(d$spectra)
```

Rows of `metadata`, `spectra`, `prediction` and `prediction_type` are aligned. `spectra` is a double-precision matrix with `wn_` wavenumber column names (cm^-1); it contains the numeric input to the saved model pipeline, not instrument acquisition files. `spectrum_id` is a new sequential public identifier. `group_id` retains the original operation-group partition and cross-validation membership. `patient_id` uses the shared P-ID namespace where an existing mapping is available; one early group has no established P-ID mapping and is left missing. No patient identity is inferred from that absence. All late bags retain their existing P IDs. `fold` is present for early OOF observations only. `y` is the saved binary target. CWM is the display label for the source NAWM category.

Preprocessing and model fitting used early CE/CWM observations only. Early IM/ED OOF bags supplied the threshold. Late observations were used for prediction only. The original fold memberships and post-fold RNG state are retained; no replacement grouping or split is generated here.

`methods.R` defines the saved computational methods adapted to the public field names. It does not execute training when sourced. It retains the original balanced-accuracy threshold rule and tie ordering, but the supplied threshold is fixed. Training functions were not run during packaging. For the published figure use the supplied predictions and threshold together, or the verified saved-model inference above.
