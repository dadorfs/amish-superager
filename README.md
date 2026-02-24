# Amish Superager Study

Code repository for the Amish Superager Study. This repository contains
analysis scripts associated with a published analysis comparing superagers
(SA), cognitively unimpaired individuals (CU), and individuals with
Alzheimer's disease (AD) in an Amish cohort.

Data are not included in this repository due to privacy restrictions.
File paths in all scripts are generic placeholders and will need to be
updated to match your local directory structure before running.

## Repository structure

```
.
└── association/    # Single-variant association analysis within linkage regions
```

## Phenotype groups

| Abbreviation | Description |
|---|---|
| SA | Superager |
| CU | Cognitively unimpaired (age-matched, non-superager) |
| AD | Alzheimer's disease dementia |

## Dependencies

All analysis scripts are written in R. Key packages:

- [GENESIS](https://bioconductor.org/packages/GENESIS/) — mixed model association testing
- [SeqArray](https://bioconductor.org/packages/SeqArray/) / [SeqVarTools](https://bioconductor.org/packages/SeqVarTools/) — GDS file handling
- [SNPRelate](https://bioconductor.org/packages/SNPRelate/) — kinship and PCA
- [tidyverse](https://www.tidyverse.org/) — data manipulation
- [data.table](https://rdatatable.gitlab.io/data.table/) — fast I/O
- [here](https://here.r-lib.org/) — project-relative file paths
