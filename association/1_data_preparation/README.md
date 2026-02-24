# 1_data_preparation

Scripts to prepare genetic and phenotypic data for association analysis.

## Overview

This directory covers all preprocessing steps needed before running association tests: converting imputed genotype data to GDS format, extracting directly genotyped variants for kinship/PCA, computing ancestry principal components and kinship coefficients, and building phenotype files with binary case-control outcomes.

## Workflow

Scripts are numbered in execution order. Phenotype preparation (`4_prepare_pheno.R`) is independent of the genotype pipeline and can be run at any point.

```
1_imputed_vcf_to_gds.R          Convert imputed VCFs → GDS
        ↓
2_create_genotyped_gds.R        Extract typed variants → merged GDS
        ↓
3_kinship_and_pcs.R             Kinship (KING) + PCA (PC-AiR, PC-Relate)

4_prepare_pheno.R               Build case-control phenotype files (independent)
```

## Scripts

| Script | Purpose | Key Output |
|--------|---------|------------|
| `1_imputed_vcf_to_gds.R` | Convert per-chromosome imputed VCF files to SeqArray GDS format, importing dosage scores (DS) | `chr*.r205.gds` |
| `2_create_genotyped_gds.R` | Subset each imputed GDS to directly genotyped variants (TYPED flag), then merge chromosomes into a single genome-wide file | `merged.maf01.r206.typed.vars.gds` |
| `3_kinship_and_pcs.R` | LD-prune genotyped variants; estimate KING-robust kinship; run PC-AiR for ancestry PCs and PC-Relate for ancestry-adjusted kinship | `pcs.RDS`, `pcrel.RDS` |
| `4_prepare_pheno.R` | Create binary case-control phenotype files for three comparisons (see below) | `sa_vs_sa_cu80plus.csv`, `sa_vs_ad.csv`, `sa_cu80plus_vs_ad.csv` |

## Phenotype Comparisons

| File | Cases (STATUS_BIN = 1) | Controls (STATUS_BIN = 0) |
|------|------------------------|---------------------------|
| `sa_vs_sa_cu80plus.csv` | Superagers (SA_OVERALL = 1, CDX = CU) | Cognitively unimpaired aged 80+ (SA_OVERALL = 0, CDX = CU, age ≥ 80) |
| `sa_vs_ad.csv` | Superagers | AD dementia (CDX = AD Dementia) |
| `sa_cu80plus_vs_ad.csv` | Superagers + CU aged 80+ | AD dementia |

## Dependencies

All scripts use the `here` package for portable path management. Set the `here` root to your project directory before running.

**R packages:** `tidyverse`, `here`, `SeqArray`, `SeqVarTools`, `SNPRelate`, `GENESIS`, `BiocParallel`
