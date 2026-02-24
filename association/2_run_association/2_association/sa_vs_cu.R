## -------------------------------------------------
##
## Script name: Single-variant association: superagers vs. cognitively unimpaired
##
## Purpose: Run single-variant association tests comparing superagers
##   (SA) to age-matched cognitively unimpaired non-superagers (CU) within
##   1-LOD linkage intervals. A logistic mixed model is fit using age, sex,
##   and top PCs as covariates, with a kinship matrix from PC-Relate to
##   account for family structure.
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
library(SeqVarTools)
library(SNPRelate)
library(Biobase)
library(GENESIS)
library(BiocParallel)
library(batchtools)
library(here)

## -------------------------------------------------
## setup multicore backend

bpparam <- MulticoreParam(workers = 10)

## -------------------------------------------------
## load data

# pheno
pheno <- read_csv(here("data", "pheno", "sa_vs_cu.csv"))
pheno$sample.id <- pheno$AGDBID

# principal components from PC-AiR
pcs <- read_rds(here("data", "pca", "pcs.RDS")) %>%
  rownames_to_column(var = "sample.id")

# PC-Relate kinship object
pcrel <- read_rds(here("data", "pca", "pcrel.RDS"))

# 1-LOD regions of interest from linkage analysis
roi <- read_csv(here("data", "genomic_ranges", "1lod_regions_regions_overall.csv"))


## -------------------------------------------------
## process

# join PCs to pheno
pheno_pcs <- pheno %>%
  mutate(sample.id = as.character(sample.id)) %>%
  left_join(pcs, by = "sample.id")

# wrap in AnnotatedDataFrame for GENESIS compatibility
annot <- AnnotatedDataFrame(pheno_pcs)

# convert PC-Relate output to a kinship matrix (scaleKin = 2 for standard kinship scale)
kinship_mat <- GENESIS::pcrelateToMatrix(pcrel, scaleKin = 2)


## -------------------------------------------------
## null model

# logistic mixed model: age, sex, and top 2 PCs as covariates;
# kinship matrix accounts for family structure
nullmod <- fitNullModel(
  annot,
  outcome = "STATUS_BIN",
  family = "binomial",
  covars = c("AGE_AT_EXAM",
             "SEX",
             paste0("PC", 1:2)),
  cov.mat = kinship_mat)

## -------------------------------------------------
## single variant testing

# run.assoc: runs single-variant association for a single GDS file
#   gds_file:   path to per-chromosome imputed GDS
#   out_prefix: prefix for the output filename
#   subdir:     subdirectory under assoc_output/ for results
#   window:     if TRUE, restrict testing to the 1-LOD region of interest for that chromosome
run.assoc <- function(gds_file, out_prefix, subdir, window = FALSE) {

  chrom <- basename(gds_file) %>%
    str_replace(".*(chr\\d*).*$", "\\1")

  # load gds
  gds <- seqOpen(gds_file)

  if (window) {

    chr_num <- as.numeric(sub("chr", "", chrom))

    lower <- filter(roi, chr == chr_num) %>%
      pull(start)

    upper <- filter(roi, chr == chr_num) %>%
      pull(end)

    seqSetFilterChrom(gds, include = chr_num, from.bp = lower, to.bp = upper)

  }

  # ensure output directory exists
  out_dir <- here("data", "assoc_output",
                  subdir)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # create SeqVarData and iterator objects
  seqData <- SeqVarData(gds)
  iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)

  # run single-variant association test
  assoc <- assocTestSingle(iterator,
                           nullmod,
                           imputed = TRUE,
                           verbose = FALSE,
                           BPPARAM = bpparam)

  # reset filter to access full variant annotation for ID lookup
  seqResetFilter(gds)

  # retrieve variant IDs (e.g. rsIDs) and merge with association results
  join_info <- seqGetData(gds, c("variant.id", "annotation/id")) %>%
    as_tibble() %>%
    rename_with(., ~ tolower(str_replace(.x, "^.*/(.*)$", "\\1")))

  assoc_annot <- assoc %>%
    left_join(join_info) %>%
    relocate(id, .after = variant.id)

  seqClose(gds)

  # export
  out_fn <- paste(chrom, out_prefix, "csv", sep = ".")

  write_csv(assoc_annot, file.path(out_dir, out_fn))

}

# list and sort imputed GDS files by chromosome
gds_files <- list.files(here("data", "imputed_gds"),
                        pattern = "chr\\d*\\.r2\\.filtered",
                        full.names = TRUE)

gds_chroms <- as.numeric(sub(pattern = ".*chr(\\d*).*$", "\\1", gds_files))

gds_files <- gds_files[order(gds_chroms)]

# run association for each chromosome in the 1-LOD regions of interest
invisible(lapply(gds_files[roi$chr],
                 function(x) run.assoc(x, out_prefix = "assoc.sa.vs.cu", subdir = "sa_vs_cu", window = TRUE)))
