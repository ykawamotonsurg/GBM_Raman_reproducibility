# GBM_Raman_reproducibility

This repository provides deidentified processed datasets and figure-specific R scripts for reproducing selected manuscript figures from our glioblastoma Raman spectroscopy study.

The current public package is aligned to the Neuro-Oncology submission-ready figure layout and includes:
- Fig. 1B-D
- Fig. 2A-D
- Fig. 3B-C

These materials are intended for **figure-level reproduction** and **manuscript transparency**. The repository does **not** include raw patient identifiers, surgery-date identifiers, raw spectra archives, raw microscopy files, raw Xenium objects, spatial coordinates, or local file-system paths.

## Repository structure

```text
GBM_Raman_reproducibility/
  README.md
  fig1_public_data/
  fig2_public_data/
  fig3_public_data/
```

Each figure folder contains processed CSV files and figure-specific R scripts used to reproduce the corresponding manuscript panels. Outputs are written to an `outputs/` subfolder inside each figure folder when the scripts are run.

## Scope

This repository is designed to support:
- reproduction of selected manuscript figure panels from processed public tables
- transparency of downstream figure-level summary data
- reuse of figure-specific plotting scripts

It is not intended to provide raw source data or internal exploratory analysis files.

## Deidentification

All public tables were prepared in deidentified form.

- Original case identifiers were replaced with public IDs such as `Case01`, `Case02`, ...
- No surgery dates, original operation identifiers, raw coordinates, or internal crosswalk tables are included
- Column names were preserved where needed to maintain compatibility with the public figure scripts

## Requirements

The figure scripts were prepared in R. Commonly used packages include:
- `readr`
- `dplyr`
- `ggplot2`
- `tidyr`
- `forcats`

Additional packages may be required depending on the script.

## How to use

1. Download or clone this repository.
2. Run the corresponding R script from the repository root or from the relevant figure folder.
3. The output figures are written to the figure folder's `outputs/` directory.

Example:

```bash
Rscript fig1_public_data/fig1b_make_public.R
Rscript fig1_public_data/fig1c_make_public.R
Rscript fig1_public_data/fig1d_make_public.R
Rscript fig2_public_data/fig2a_make_public.R
Rscript fig2_public_data/fig2b_make_public.R
Rscript fig2_public_data/fig2c_make_public.R
Rscript fig2_public_data/fig2d_make_public.R
Rscript fig3_public_data/fig3b_program_scores_ir_msr1_state_make_public.R
Rscript fig3_public_data/fig3c_msr1pos_programs_ir_percase_make_public.R
```

## Notes

This package intentionally contains only main-figure reproduction scripts and processed public tables. Supplementary diagnostic or radiographic-context analyses are not included in this main-figure code package.
