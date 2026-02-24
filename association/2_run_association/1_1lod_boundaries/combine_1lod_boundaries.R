## -------------------------------------------------
##
## Script name: Identify peak 1-LOD boundaries across dominant and recessive analyses
##
## Purpose: Combine 1-LOD boundary intervals from the dominant and
##   recessive parametric linkage analyses into a single set of genomic regions.
##   For chromosomes where a significant HLOD was observed under both models,
##   the boundary from the model with the higher peak HLOD is retained.
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

rec_1lod <- read_csv(here("data", "genomic_ranges", "pl_rec_1lod_regions_of_interest.csv"))

rec_out <- read_tsv(here("data", "linkage", "pl_rec.out"))

dom_1lod <- read_csv(here("data", "genomic_ranges", "pl_dom_1lod_regions_of_interest.csv"))

dom_out <- read_tsv(here("data", "linkage", "pl_dom.out"))

## -------------------------------------------------
## process

# tag each boundary table by model before combining
rec_1lod$run <- "rec"
dom_1lod$run <- "dom"

lod_down_joined <- rbind(rec_1lod, dom_1lod)

# extract peak HLOD per chromosome for each model
dom_max_lod <- dom_out %>%
  slice_max(by = CHR, n = 1, order_by = HLOD) %>%
  select(CHR, HLOD) %>%
  mutate(run = "dom")

rec_max_lod <- rec_out %>%
  slice_max(by = CHR, n = 1, order_by = HLOD) %>%
  select(CHR, HLOD) %>%
  mutate(run = "rec")

max_lod_all <- bind_rows(dom_max_lod, rec_max_lod)

# join peak HLOD back to combined boundary table
lod_down_joined_maxlod <- lod_down_joined %>%
  left_join(max_lod_all, join_by(chr == CHR, run))

# for each chromosome, retain the boundary from the model with the higher HLOD
lod_down_filter <- lod_down_joined_maxlod %>%
  slice_max(by = chr, n = 1, order_by = HLOD) %>%
  arrange(chr)

## -------------------------------------------------
## export

out <- lod_down_filter %>%
  select(chr, start, end)

write_csv(out,
          here("data", "genomic_ranges", "1lod_regions_overall.csv"))


