# 2_run_association

Scripts for targeted single-variant association analysis within linkage-derived genomic intervals.

## Overview

Rather than testing the full genome, association testing is restricted to regions identified by parametric linkage analysis (1-LOD intervals around significant HLOD peaks).

Association tests compare superagers (SA) to two phenotypic groups — individuals with Alzheimer's disease (AD) and cognitively unimpaired non-superagers (CU) — each with and without APOE e4 adjustment. Multiple testing correction uses the SimpleM method, which accounts for linkage disequilibrium among tested variants.

## Workflow

```
Step 1 → Define 1-LOD boundary intervals from linkage peaks
Step 2 → Run single-variant association tests within those intervals
Step 3 → Compute multiple testing corrections (SimpleM)
```

## Directory Structure

```
2_run_association/
│
├── 1_1lod_boundaries/
│   ├── calculate_dominant_1lod.R      # 1-LOD intervals from dominant HLOD peaks (chr 1, 2, 7, 20)
│   ├── calculate_recessive_1lod.R     # 1-LOD intervals from recessive HLOD peaks (chr 2, 16, 20)
│   └── combine_1lod_boundaries.R      # Merge dominant/recessive intervals; retain higher-HLOD boundary per chromosome
│
├── 2_association/
│   ├── sa_vs_ad.R                     # Superagers vs. Alzheimer's disease
│   ├── sa_vs_ad_apoe_adjusted.R       # Superagers vs. AD, adjusted for APOE e4 carrier status
│   ├── sa_vs_cu.R                     # Superagers vs. cognitively unimpaired
│   └── sa_vs_cu_apoe_adjusted.R       # Superagers vs. CU, adjusted for APOE e4 carrier status
│
└── 3_multiple_testing/
    ├── create_simplem_matrices.R      # Extract per-chromosome dosage matrices within 1-LOD windows
    └── run_simplem.R                  # Estimate effective number of tests (Meff) and adjusted alpha per region
```

## Key Inputs / Outputs

| Step | Key Inputs | Key Outputs |
|------|-----------|-------------|
| 1a – dominant boundaries | `pl_dom.out` | `pl_dom_1lod_regions_of_interest.csv` |
| 1b – recessive boundaries | `pl_rec.out` | `pl_rec_1lod_regions_of_interest.csv` |
| 1c – combined boundaries | above two CSVs | `1lod_regions_overall.csv` |
| 2 – association | `1lod_regions_overall.csv`, imputed GDS, kinship matrix, PCs | per-chromosome association result CSVs |
| 3a – SimpleM matrices | imputed GDS, `1lod_regions_overall.csv` | `chr*_1lod_simplem_mtx.csv` |
| 3b – SimpleM correction | above matrices | `simplem_tbl_1lod.csv` |

## Models

All association tests use a **logistic mixed model** (GENESIS `fitNullModel`) with:
- **Outcome:** binary SA status
- **Covariates:** age at exam, sex, top 2 PCs (APOE e4 carrier status added in adjusted models)
- **Random effect:** kinship matrix from PC-Relate to account for family structure

## Dependencies

All scripts use the `here` package for portable path management. Set the `here` root to your project directory before running.

**R packages:** `tidyverse`, `data.table`, `here`, `SeqArray`, `SeqVarTools`, `SNPRelate`, `Biobase`, `GENESIS`, `BiocParallel`, `batchtools`
