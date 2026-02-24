# Amish Superager Association Analysis

Code for the genetic association analysis component of the Amish Superager Study. The analysis identifies genetic variants associated with the superager phenotype by restricting single-variant tests to genomic regions implicated by prior parametric linkage analysis.

## Study Background

Superagers are individuals aged 80 and older who maintain cognitive performance at or above the level of adults decades younger. This analysis tests for genetic variants associated with superager status in the Amish population, comparing superagers to individuals with Alzheimer's disease (AD) and to cognitively unimpaired (CU) non-superagers.

Family structure is accounted for using a kinship matrix derived from PC-Relate. Association testing is restricted to 1-LOD intervals around significant HLOD peaks from parametric linkage analysis, reducing the multiple testing burden while focusing power on linkage-supported regions.

## Workflow

```
1_data_preparation/     Prepare genotype and phenotype data for analysis
        ↓
2_run_association/      Define 1-LOD boundaries and run single-variant association tests
```

## Directory Structure

```
association/
├── 1_data_preparation/
│   ├── 1_imputed_vcf_to_gds.R          Convert imputed VCFs to GDS format
│   ├── 2_create_genotyped_gds.R         Extract typed variants; merge into genome-wide GDS
│   ├── 3_kinship_and_pcs.R              KING kinship; PC-AiR ancestry PCs; PC-Relate kinship
│   └── 4_prepare_pheno.R                Build binary case-control phenotype files
│
└── 2_run_association/
    ├── 1_1lod_boundaries/
    │   ├── calculate_dominant_1lod.R    1-LOD intervals from dominant HLOD peaks
    │   ├── calculate_recessive_1lod.R   1-LOD intervals from recessive HLOD peaks
    │   └── combine_1lod_boundaries.R    Merge intervals; retain higher-HLOD boundary per chromosome
    ├── 2_association/
    │   ├── sa_vs_ad.R                   Superagers vs. AD
    │   ├── sa_vs_ad_apoe_adjusted.R     Superagers vs. AD, APOE e4-adjusted
    │   ├── sa_vs_cu.R                   Superagers vs. cognitively unimpaired
    │   └── sa_vs_cu_apoe_adjusted.R     Superagers vs. CU, APOE e4-adjusted
    └── 3_multiple_testing/
        ├── create_simplem_matrices.R    Extract per-chromosome dosage matrices
        └── run_simplem.R                Compute effective number of tests (SimpleM)
```

See the README in each subdirectory for details on inputs, outputs, and dependencies.

## Association Models

All tests use a **logistic mixed model** (GENESIS `fitNullModel`) with:
- **Outcome:** binary superager status
- **Covariates:** age at exam, sex, top 2 ancestry PCs
- **APOE-adjusted models:** additionally include APOE e4 carrier status
- **Random effect:** PC-Relate kinship matrix to account for family structure

## Dependencies

All scripts use the `here` package for path management. File paths in each script are generic placeholders.

**R packages:** `tidyverse`, `data.table`, `here`, `SeqArray`, `SeqVarTools`, `SNPRelate`, `Biobase`, `GENESIS`, `BiocParallel`, `batchtools`
