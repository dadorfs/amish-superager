## -------------------------------------------------
##
## Script name: Create genotype matrices for SimpleM multiple testing correction
##
## Purpose: For each chromosome in the 1-LOD regions of interest,
##   extract the imputed dosage matrix from GDS files. These per-chromosome
##   matrices are used as input to run_simplem.R to estimate the effective
##   number of independent tests (Meff) for multiple testing correction.
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
library(data.table)
library(here)
library(SeqArray)
library(SeqVarTools)

## -------------------------------------------------
## load data

# pheno
pheno <- read_csv(here("data", "pheno", "sa_vs_ad.csv"))
pheno$sample.id <- pheno$CGI

# osa 1-lod intervals
lod_down <- read_csv(here("data", "genomic_ranges", "1lod_regions_regions_overall.csv"))


## -------------------------------------------------
## process data 

## -------------------------------------------------
## create simplem function 

# function to create genotype matrix for osa windows
create.gt.matrix <- function(chrom){

    pattern <- sprintf("chr%d.r2.filtered.gds", chrom)
    gds_dir <- here("data", "imputed_gds")

    gds_files <- list.files(gds_dir, pattern = pattern, full.names = TRUE)
    gds_file <- gds_files[1]

    lower <- filter(lod_down, chr == chrom) %>%
        pull(start) 

    upper <- filter(lod_down, chr == chrom) %>%
        pull(end) 

    gds <- seqOpen(gds_file)

    seqSetFilter(gds, sample.id = pheno$AGDBID)
    seqSetFilterChrom(gds, include = chrom, from.bp = lower, to.bp = upper)
    seqSetFilterCond(gds, maf = 0.01)

    dosage_mtx <- t(seqGetData(gds, "$dosage_alt"))

    out_dir <- here("data", "simplem")


    out_fn <- sprintf("chr%d_1lod_simplem_mtx.csv", chrom)

    fwrite(dosage_mtx, file.path(out_dir, out_fn))

    seqClose(gds)

}

## -------------------------------------------------
## create simplem matrices 
    
lapply(lod_down$chr, function(x) create.gt.matrix(x))













