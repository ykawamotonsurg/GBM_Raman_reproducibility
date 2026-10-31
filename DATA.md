# Data guide

| Result | Analysis unit and input | Code / location |
|---|---|---|
| Fig. 1 | Main 57-bag scores, saved ROC points/statistics, mean spectra and patient-by-region peak means | `main/code/plot_saved.py` |
| Fig. 2 | 1,499 included cells; fixed UMAP, two peak measurements, scores and 225 high labels; 38 case-region summaries | `shared/data/cells.csv`; fixed threshold and standardization files; `main/code/summarize_fixed_inputs.py` |
| Fig. 3B / Table S6 | 76 same-cell records, including 55 evaluable and 21 missing MSR1 measurements | `main/data/samecell_labels.csv`; `main/code/summarize_fixed_inputs.py` |
| Fig. 3C | Measurement-unit means, retaining unequal group sizes; THP1 control n=3 and loaded n=7 | `main/data/culture_unit_means.csv`; `main/code/summarize_fixed_inputs.py` |
| Fig. 3D / Table S7 | Six paired THP1 units on three experimental days, with two pairs per day | `main/data/Fucoidan_Pairs.csv`; saved statistics in `Table_S7A.csv` and `Table_S7B.csv` |
| Fig. 4B / Figs. S4–S6 | Case/niche/state/program scores, primary three-case values and four-case sensitivity; program gene membership | `shared/data/spatial_program_scores.csv`, `spatial_gene_sets.csv`, `supplement/data/spatial_program_scores_all4.csv`; supplement plotting code |
| Fig. S2B | Fixed temporal development/evaluation spectra, folds, predictions, preprocessing and model | `supplement/analysis/s2b/` |
| Fig. S3 / Table S5 | Fixed per-cell inputs and 37 matched tissue–cell pairs | `shared/data/cells.csv`, `supplement/data/tissue_cell_bridge.csv`; `supplement/analysis/statistics/fixed_input_statistics.R` |
| Fig. S7 / Table S9 | 18 MRI–cell cases; shared two-rater measured volumes and fixed CE cell counts | `supplement/data/MRI_CaseData.csv`; `supplement/analysis/statistics/mri_agreement.R`; `fixed_input_statistics.R` |

P IDs are reused across established patient-linked tables. Same-cell IDs belong to a separate experiment and must not be joined to the UMAP cell IDs. Culture `unit_id` preserves pairing and `experiment_id` preserves the day block. Spatial GBM-A/B/C/D identifiers belong to a separate namespace; no cross-modality patient join is implied. Empty entries remain missing, not zero or negative.

The article's **Data S1** provides the additional tabulated results: `Table_S1`–`Table_S9B`, `Singlecell_Trend`, `Singlecell_Adjusted`, `Influence_FullPrecision`, `QC_MSR1`, `Xenium_StateMeans` and `Xenium_D_niches`. These are source-data supplement sheet names. `QC_MSR1` describes the culture experiment and is not the source of the same-cell labels. The repository does not duplicate every supplementary table.

The plotting commands use the supplied fixed results. The optional statistical code documents the saved trend/association/MRI methods from case-level inputs; it does not recompute UMAP or cell labels. MRI association and rater-agreement statistics have been regenerated from the shared measured-volume input. The non-MRI statistical scopes have not been rerun. The S2B model-inference path has been checked after export.

`MRI_CaseData.csv` contains total and necrotic volumes for both raters in mm³, alongside the fixed CE high/total cell counts. It is the common measurement input for both MRI panels. Its `enhance_fraction` equals `(rater1_total_volume_mm3 - rater1_necrotic_volume_mm3) / rater1_total_volume_mm3`.

For Fig. S7A, `mri_agreement.R` derives enhancing volumes and fractions for both raters, then regenerates `mri_agreement.csv` and the paired-volume mean/difference tables `mri_tumor_volume.csv` and `mri_necrotic_volume.csv`. These saved tables are the plotter's inputs. For Fig. S7B, the plotter uses `MRI_CaseData.csv` and reads n, rho and P from `mri_spearman.csv`. The MRI-only `fixed_input_statistics.R` command derives the association fraction from the same rater-1 volumes, verifies the stored fraction, and regenerates the results underlying Table S9B and its supporting analyses. Plotting performs no statistical tests.

Instrument-level tissue processing and the main classifier training pipeline, the original UMAP fitting inputs, microscopy-to-cell mapping, and the initial spatial count/segmentation-to-score workflow are not reconstructed in this package. The spatial gene-difference table contains saved summaries; the original case-to-gene aggregation code is not supplied. The original paired-test and tissue-peak inference implementations are not included. These limits should not be described as full raw-data analysis reproduction.
