library(tidyverse)

#DIVERSITY_IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/diversity_metrics/CD4_Mem-Sorted_diversity.rds"
DIVERSITY_IN_PATH <- snakemake@input[[1]]

DIVERSITY_OUT_PATH <- snakemake@output[[1]]

KEEP_TIMEPOINTS <- snakemake@params[["keep_timepoints"]]

diversity_dat <- readRDS(DIVERSITY_IN_PATH)

diversity_dat <- diversity_dat %>%
        filter(Timepoint %in% KEEP_TIMEPOINTS)

saveRDS(diversity_dat, DIVERSITY_OUT_PATH)
