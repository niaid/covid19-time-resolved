library(tidyverse)

PBULKLIST_IN_PATH <- snakemake@input[[1]]

PBULKLIST_OUT_PATH <- snakemake@output[[1]]

pbulklist <- readRDS(PBULKLIST_IN_PATH)

pbulklist <- lapply(pbulklist, function(dge){
  meta <- dge$samples
  #Add onset group
  meta <- meta %>%
          mutate(onset_group = NA) %>%
          mutate(onset_group = replace(onset_group, days_since_onset < 14, "early")) %>%
          mutate(onset_group = replace(onset_group, days_since_onset >= 14 & days_since_onset <= 21, "mid")) %>%
          mutate(onset_group = replace(onset_group, days_since_onset > 21, "late"))

  meta %>% select(Donor, Timepoint, Batch, onset_group, days_since_onset) %>% distinct

  #add severity.outcome column that has HC for healthy controls
  meta <- meta %>%
          mutate(severity.outcome2 = as.character(severity.outcome)) %>%
          mutate(severity.outcome2 = replace(severity.outcome2, Class =="HC", "HC")) %>%
          mutate(severity.outcome2 = relevel(factor(severity.outcome2), "HC"))

  dge$samples <- meta
  return(dge)
})

saveRDS(pbulklist, PBULKLIST_OUT_PATH)

