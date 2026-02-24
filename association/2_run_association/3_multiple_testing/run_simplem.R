## -------------------------------------------------
##
## Script name: Run SimpleM multiple testing correction
##
## Purpose: Apply the SimpleM method (Gao et al. 2008, 2009) to
##   estimate the effective number of independent tests (Meff) within each
##   1-LOD region of interest. Uses per-chromosome genotype matrices produced
##   by create_simplem_matrices.R. Outputs adjusted significance thresholds
##   (alpha_adj = 0.05 / Meff) per region.
##
## -------------------------------------------------
##
## Notes:
##   File paths throughout this script are generic placeholders.
##   Update paths to match your local directory structure before running.
##
##   SimpleM citation:
##   Gao X, Starmer J, Martin ER (2008). A Multiple Testing Correction Method
##   for Genetic Association Studies Using Correlated Single Nucleotide
##   Polymorphisms. Genetic Epidemiology 32:361-369.
##
##   Gao X, Becker LC, Becker DM, Starmer J, Province MA (2009). Avoiding the
##   high Bonferroni penalty in genome-wide association studies. Genetic
##   Epidemiology (Epub ahead of print).
##
## -------------------------------------------------

library(tidyverse)
library(data.table)
library(here)

#============================================================================
# Meff through the PCA approach
# use a part of the eigen values according to how much percent they contribute
# to the total variation 

Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if(myEigenSum <= myCut){
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    }
    else{
      break
    }
  }	
  return(index_Eigen)
}

#============================================================================
# infer the cutoff => Meff

inferCutoff <- function(dt_My){
  CLD <- cor(dt_My)
  eigen_My <- eigen(CLD)
  
  # PCA approach
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}

#============================================================================

PCA_cutoff <- 0.995

#============================================================================
# fix length, simpleM

run.simplem <- function(chrom){

    pattern <- sprintf("chr%d_1lod_simplem_mtx.csv", chrom)

    simplem_dir <- here("data", "simplem")

    simplem_files <- list.files(simplem_dir, pattern = pattern, full.names = TRUE)
    simplem_file <- simplem_files[1]

    mySNP_nonmissing <- read.csv(simplem_file, colClasses="integer")		

    numLoci <- length(mySNP_nonmissing[, 1])

    simpleMeff <- NULL
    fixLength <- 133 
    i <- 1
    myStart <- 1
    myStop <- 1
    while(myStop < numLoci){
        myDiff <- numLoci - myStop 
        if(myDiff <= fixLength) break

        myStop <- myStart + i*fixLength - 1
        snpInBlk <- t(mySNP_nonmissing[myStart:myStop, ])
        MeffBlk <- inferCutoff(snpInBlk)
        simpleMeff <- c(simpleMeff, MeffBlk)
        myStart <- myStop+1
    }
    snpInBlk <- t(mySNP_nonmissing[myStart:numLoci, ])
    MeffBlk <- inferCutoff(snpInBlk)
    simpleMeff <- c(simpleMeff, MeffBlk)

    # cat("Total number of SNPs is: ", numLoci, "\n")
    # cat("Inferred Meff is: ", sum(simpleMeff), "\n")

    return(list(numLoci = numLoci, simpleMeff = sum(simpleMeff)))


}

#============================================================================
# run simplem
# osa 1-lod intervals
lod_down <- read_csv(here("data", "genomic_ranges", "1lod_regions_regions_overall.csv"))


results <- lapply(lod_down$chr, run.simplem)

simplem_tbl <- tibble(
                      chromosome = lod_down$chr,
                      numLoci = sapply(results, function(x) x$numLoci),
                      simpleMeff = sapply(results, function(x) x$simpleMeff))

simplem_tbl <- simplem_tbl %>%
    mutate(alpha_adj = 0.05 / simpleMeff)

#============================================================================
# join boundaries
simplem_tbl_bnds <- simplem_tbl %>%
    left_join(lod_down, join_by(chromosome == chr)) %>%
    relocate(start, end, .after = chromosome)

#============================================================================
# export 
write_csv(simplem_tbl_bnds,
          here("data", "simplem", "simplem_tbl_1lod.csv"))

