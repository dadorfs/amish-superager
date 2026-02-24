## -------------------------------------------------
##
## Script name: Compute 1-LOD boundaries for parametric dominant linkage
##
## Purpose: For each chromosome with a significant dominant HLOD
##   peak, identify the genomic region within 1 LOD unit of the peak. These
##   intervals define the regions used for downstream single-variant
##   association analysis.
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

## -------------------------------------------------
## load data

pl_out <- fread(here("data", "linkage", "pl_dom.out"))

## -------------------------------------------------
## process

# convert LABEL to bp position
pl_out[, LABEL := as.numeric(sub(".*:(\\d+)-.*", "\\1", LABEL))]

# add maximum hlod and peak position columns
pl_out[, MAX_HLOD := max(HLOD), by = CHR]
pl_out[, PEAK_POS := LABEL[which.max(HLOD)], by = CHR]

pl_1lod_filter <- data.table()

# chromosomes with a significant dominant HLOD peak
for (chr in c(1, 2, 7, 20)) {

  # filter to current chromosome
  chr_data <- pl_out[CHR == chr]

  # record peak position and maximum HLOD for this chromosome
  peak_pos <- chr_data$PEAK_POS[1]
  max_hlod <- chr_data$MAX_HLOD[1]

  # split data into regions before and after the peak
  before_peak <- chr_data[LABEL <= peak_pos]
  after_peak  <- chr_data[LABEL >= peak_pos]

  # walk backward from the peak to find where HLOD drops by 1 unit
  start_pos <- peak_pos
  for (i in seq(nrow(before_peak), 1, -1)) {
    if ((max_hlod - before_peak$HLOD[i]) > 1) {
      start_pos <- before_peak$LABEL[i]
      break
    }
  }

  # walk forward from the peak to find where HLOD drops by 1 unit
  end_pos <- peak_pos
  for (i in seq(1, nrow(after_peak))) {
    if ((max_hlod - after_peak$HLOD[i]) > 1) {
      end_pos <- after_peak$LABEL[i]
      break
    }
  }

  # extract all markers within the 1-LOD interval
  contiguous_data <- chr_data[LABEL >= start_pos & LABEL <= end_pos]

  pl_1lod_filter <- rbind(pl_1lod_filter, contiguous_data, fill = TRUE)

}

pl_1lod_boundaries <- pl_1lod_filter[, .(start = min(LABEL), end = max(LABEL)), by = CHR]

setnames(pl_1lod_boundaries, tolower(names(pl_1lod_boundaries)))


## -------------------------------------------------
## export

fwrite(pl_1lod_boundaries, 
       here("data", "genomic_ranges", "pl_dom_1lod_regions_of_interest.csv"))



