## -------------------------------------------------
##
## Script name: Prepare phenotype files for association analysis
##
## Purpose: Build binary outcome variables (STATUS_BIN) for three
##   case-control comparisons: superagers vs. cognitively unimpaired 80+,
##   superagers vs. AD dementia, and superagers+CU80+ vs. AD dementia.
##   Only genotyped individuals are retained.
##
## -------------------------------------------------
##
## Notes:
##   File paths throughout this script are generic placeholders.
##   Update paths to match your local directory structure before running.
##
## -------------------------------------------------

## load packages

library(tidyverse)
library(here)

## -------------------------------------------------
## load data

sa <- read_csv(here("data", "raw", "phenotype.csv"))

## -------------------------------------------------
## process data

# Select variables
sa_slct <- sa %>%
    select(
    AGDBID,             # Study ID
    SEX,
    AGE_AT_EXAM,
    SA_OVERALL,         # SuperAger status
    CDX_final,          # Consensus diagnosis
    E4_DOSAGE,
    GENOTYPED)

# Filter genotyped
sa_slct_fltr <- sa_slct %>%
    filter(GENOTYPED == 1)

# Assign status
final_sa_cu80plus <- sa_slct_fltr %>%
    mutate(STATUS_BIN = 
           case_when(
                     # SuperAger and unimpaired by consensus
                     SA_OVERALL == 1 & CDX_final == "CU" ~ 1,
                     # Non-SuperAger, cognitively unimpaired by consensus, and ≥ 80
                     SA_OVERALL == 0 & CDX_final == "CU" & AGE_AT_EXAM >= 80 ~ 0,
                     .default = NA_real_)
    )


final_sa_ad <- sa_slct_fltr %>%
    mutate(STATUS_BIN = 
           case_when(
                     # SuperAger and unimpaired by consensus
                     SA_OVERALL == 1 & CDX_final == "CU" ~ 1,
                     # AD by consensus
                     CDX_final == "AD Dementia" ~ 0,
                     .default = NA_real_)

    )

final_sa_cu80plus_ad <- sa_slct_fltr %>%
    mutate(STATUS_BIN = 
           case_when(
                     # SuperAger and unimpaired by consensus OR 
                     # Non-SuperAger, cognitively unimpaired by consensus, and ≥ 80
                     SA_OVERALL == 1 & CDX_final == "CU" ~ 1,
                     SA_OVERALL == 0 & CDX_final == "CU" & AGE_AT_EXAM >= 80 ~ 1,
                     # AD by consensus
                     CDX_final == "AD Dementia" ~ 0,
                     .default = NA_real_)

    )

## -------------------------------------------------
## export
out_dir <- here("data", "pheno")

write_csv(final_sa_cu80plus,
          file.path(out_dir, "sa_vs_sa_cu80plus.csv"))

write_csv(final_sa_ad,
          file.path(out_dir, "sa_vs_ad.csv"))

write_csv(final_sa_cu80plus_ad,
          file.path(out_dir, "sa_cu80plus_vs_ad.csv"))
