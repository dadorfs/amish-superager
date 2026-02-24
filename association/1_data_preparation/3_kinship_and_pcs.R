## -------------------------------------------------
##
## Script name: Calculate kinship and PCs
##
## Purpose: Estimate pairwise kinship using the KING-robust method,
##   then perform PC-AiR for ancestry estimation and PC-Relate for kinship
##   estimation that accounts for population structure.
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
library(GENESIS)
library(BiocParallel)
library(here)

## -------------------------------------------------
## setup multicore backend

bpparam <- MulticoreParam(workers = 10)

## -------------------------------------------------
## pruning and kinship estimation

# open gds connection
gds_fn <- here("data", "gtyped_gds", "merged.maf01.r206.typed.vars.gds")

gds <- seqOpen(gds_fn)

# prune SNPs for LD (r^2 < 0.16, MAF > 0.05) before kinship and PCA
set.seed(123)
snpset <- snpgdsLDpruning(gds,
                          maf = 0.05,
                          ld.threshold = 0.16,
                          num.thread = bpparam$workers)

snpset_ids <- unlist(snpset, use.names = FALSE)

# estimate pairwise kinship using the KING-robust method
king <- snpgdsIBDKING(gds,
                      snp.id = snpset_ids,
                      num.thread = bpparam$workers)

kingMat <- king$kinship

colnames(kingMat) <- rownames(kingMat) <- king$sample.id
kingMat[1:5, 1:5]

## -------------------------------------------------
## PC-AiR

# PC-AiR computes PCs robust to family structure by partitioning samples into
# related and unrelated sets; kin.thresh = 2^(-9/2) corresponds to between 3rd and 4th degree relatives
pca <- GENESIS::pcair(
  gds,
  kinobj = kingMat,
  kin.thresh = 2 ^ (-9 / 2),
  divobj = kingMat,
  div.thresh = 2 ^ (-9 / 2),
  snp.include = snpset_ids,
  num.cores = bpparam$workers
)

# extract top 10 PCs
pcs <- data.frame(pca$vectors[, 1:10])
colnames(pcs) <- paste0("PC", 1:10)

write_rds(pcs,
          here("data", "pca", "pcs.RDS"))

## -------------------------------------------------
## PC-Relate

# PC-Relate estimates kinship coefficients accounting for ancestry 
seqData <- SeqVarTools::SeqVarData(gds)
iterator <- SeqVarTools::SeqVarBlockIterator(seqData, verbose = FALSE)

pcrel <- pcrelate(iterator,
                  pcs = pca$vectors[, 1:3],
                  training.set = pca$unrels,
                  BPPARAM = bpparam)

write_rds(pcrel,
          here("data", "pca", "pcrel.RDS"))

# close gds
seqClose(gds)

