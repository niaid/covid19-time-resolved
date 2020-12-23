library(tidyverse)

eset <- readRDS("data/CITE5p/all_batches/differential_expression_cell_freq/2020_07_20/expressionsets/adjustedcelltype_pcnt_parent-eset.rds")

endpts <- readRDS("data/metadata/endpoints/2020_07_20/batch.sample.end.points.log10.RDS")

meta <- pData(eset)

meta %>% filter(Class == "COVID") %>%
        select(Timepoint,PC1, PLS1, Donor, initial.Spike)

subj_w_NA <- meta %>% filter(is.na(PC1) & Class == "COVID") %>% pull(Donor) %>% unique()

meta %>% filter(is.na(PC1) & Class == "COVID") %>% 
        select(Timepoint,PC1, Donor, days_since_onset, days_since_hospitalized)
subj_w_NA %in% endpts$subject_id
