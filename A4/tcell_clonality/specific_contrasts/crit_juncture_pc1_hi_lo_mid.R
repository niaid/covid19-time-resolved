library(tidyverse)
library(lmerTest)

IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/all_timepoints_covid_only/dat_list/CD8_Mem-Sorted_diversity.rds"

dat_list <- readRDS(IN_PATH)

dat <- dat_list$simpson

dat <- dat %>%
        mutate(onset_group = NA) %>%
        mutate(onset_group = replace(onset_group, days_since_onset < 17, "early")) %>%
        mutate(onset_group = replace(onset_group, days_since_onset >= 17 & days_since_onset <= 23, "mid")) %>%
        mutate(onset_group = replace(onset_group, days_since_onset > 23, "late"))

dat <- dat %>% mutate(group = paste(onset_group, PC1_cat))

form <- "median1000 ~ 0 + PC1_cat:onset_group + PC1_cat + onset_group + Age + (1|Batch)"
form <- "median1000 ~ 0 + group + Age + (1|Batch)"

fit <- lmer(form, data = dat)

fixef(fit)
contrast_vec <- rep(0, length(fixef(fit)))
names(contrast_vec) <- names(fixef(fit))

contrast_vec["groupmid PC1_high"] <- 1
contrast_vec["groupmid PC1_low"] <- -1

contest(fit)
