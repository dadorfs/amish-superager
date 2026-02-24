## -------------------------------------------------
##
## Script name: Extract genotyped variants from imputed GDS
##
## Purpose: Subset per-chromosome imputed GDS files to retain only
##   directly genotyped variants (TYPED flag), then merge into a single
##   genome-wide GDS file for use in kinship and PC estimation.
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
imputed_dir <- here("data", "imputed_gds")
imputed_fn <- list.files(imputed_dir, pattern = "chr\\d*.r2.filtered.gds$", full.names = TRUE)

gtyped_dir <- here("data", "gtyped_gds")


# sort GDS files by chromosome number
chroms <- sub(pattern = ".*chr(\\d*)\\.r2\\.filtered.gds*$",
              replacement = "\\1",
              imputed_fn) %>%
  as.numeric()

imputed_fn <- imputed_fn[order(chroms)]

# loop through each GDS, subsetting to directly genotyped (non-imputed) variants
output_files <- c()

for (file in imputed_fn) {

  base_fn <- basename(file) %>%
    tools::file_path_sans_ext()

  out_fn <- file.path(gtyped_dir,
                      paste0(base_fn, ".typed.vars.gds"))

  gds <- seqOpen(file)
  # TYPED flag identifies variants that were directly genotyped on the array
  typed_variants <- seqGetData(gds, "annotation/info/TYPED")

  seqSetFilter(gds, variant.sel = typed_variants)

  seqExport(gds, out_fn)

  seqClose(gds)

  output_files <- c(output_files, out_fn)

}

# merge per-chromosome GDS files into a single genome-wide GDS
seqMerge(gds.fn = output_files,
         out.fn = file.path(gtyped_dir,
                            "merged.r2.filtered.typed.vars.gds"))

