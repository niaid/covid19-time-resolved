library(tidyverse)

#DIVERSITY_IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/diversity_metrics/CD4_Mem-Sorted_diversity.rds"
DIVERSITY_IN_PATH <- snakemake@input[[1]]

#META_IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_meta.tsv"
META_IN_PATH <- snakemake@input[[2]]

DAT_LIST_OUT_PATH <- snakemake@output[[1]]

diversity_dat <- readRDS(DIVERSITY_IN_PATH)
meta <- read_tsv(META_IN_PATH)

dat <- left_join(diversity_dat, meta)

dat <- dat %>%
        mutate(donor_time_batch = paste(Donor, Timepoint, Batch, sep = "_"))

measures <- unique(dat[["measure"]])
names(measures) <- measures
dat_list <- lapply(measures, function(m){
  dat %>% filter(measure == m)
})

saveRDS(dat_list, DAT_LIST_OUT_PATH)
