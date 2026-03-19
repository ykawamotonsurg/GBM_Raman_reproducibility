# GBM_Raman_reproducibility

This repository provides deidentified processed datasets and figure-specific R scripts for reproducing selected manuscript figures from our glioblastoma Raman spectroscopy study.

The current public package includes:
- Fig. 1
- Fig. 2
- Fig. 3

These materials are intended for **figure-level reproduction** and **manuscript transparency**. The repository does **not** include raw patient identifiers, surgery-date identifiers, raw spectra archives, raw microscopy files, raw Xenium objects, spatial coordinates, or local file-system paths.

## Repository structure

```text
GBM_Raman_reproducibility/
  README.md
  fig1_public_data/
  fig2_public_data/
  fig3_public_data/
```

Each figure folder contains processed CSV files and figure-specific R scripts used to reproduce the corresponding manuscript panels.

## Scope

This repository is designed to support:
- reproduction of selected manuscript figures from processed public tables
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

Additional packages may be required depending on the script.

## How to use

1. Download or clone this repository.
2. Open the relevant figure folder.
3. Run the corresponding R script using the processed CSV files in that folder.
4. The output figures should reproduce the corresponding manuscript panels from the public processed data.

## Notes

The repository currently contains public reproduction packages for Fig. 1–3. Additional figure packages may be added in the final release if needed.
