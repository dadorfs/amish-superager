## -------------------------------------------------
##
## Script name: Convert imputed VCFs to GDS format
##
## Purpose: Convert per-chromosome imputed VCF files to SeqArray
##   GDS format, importing dosage scores (DS) for downstream association
##   analysis.
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
library(SeqArray)
library(here)

## file paths
vcf_files <- list.files(here("data", "imputed_vcf"),
                        pattern = "chr\\d*.r2.filtered.vcf.gz$",
                        full.names = TRUE)

# extract chromosome numbers and sort files in order (chr1, chr2, ..., chr22)
chrom_numbers <- as.numeric(gsub(".*chr([0-9]*).*", "\\1", vcf_files))
vcf_files_sorted <- vcf_files[order(chrom_numbers)]

## convert to GDS
# DS (dosage) field is imported; INFO fields are skipped (info.import = character(0))
lapply(vcf_files_sorted[1:22], function(f) {
  gds_file <- file.path(
    here(), "data", "imputed_gds",
    sub(".vcf.gz", ".gds", basename(f))
  )
  seqVCF2GDS(
    vcf.fn = f,
    out.fn = gds_file,
    fmt.import = "DS",
    info.import = character(0),
    parallel = 4L
  )
})
