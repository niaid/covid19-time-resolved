library(tidyverse)

IN_PATH <- "/hpcdata/sg/sg_data/users/martinsaj/COVID19/externalData/SchulteSchrepping/seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"

OUT_PATH <- "data/externalData/SchulteSchrepping/seurat_COVID19_PBMC_cohort1_add_timepoint_2020-10-30.rds"
dir.create(dirname(OUT_PATH), recursive = TRUE)

LABEL_TRANSFER_PATH <- "/hpcdata/sg/sg_data/users/martinsaj/COVID19/dataupload/SeuratObjects/predictedLabels.Schulte.Rds"

label_obj <- readRDS(LABEL_TRANSFER_PATH)

obj <- readRDS(IN_PATH)

stopifnot(identical(rownames(label_obj), colnames(obj)))
head(rownames(label_obj) == colnames(obj))
head(cbind(rownames(label_obj), colnames(obj)))

obj$transferred_labels <- label_obj$predicted.id

meta <- obj@meta.data

keep_cols <- c(
  "orig.ident",
  "donor", "sampleID", "origID", "days_after_onset", 
  "DaysAfterSample0", "experiment", "age", "sex", 
  "group_per_sample", "who_per_sample", "disease_stage", 
  "outcome", "comorbidities")

keep_cols %in% colnames(meta)

meta_sub <- meta %>% select(keep_cols) %>% distinct()

meta_sub <- meta_sub %>%
        #mutate(Timepoint = "") %>%
        group_by(donor) %>%
        mutate(Timepoint = paste0("T", rank(days_after_onset)-1)) %>%
        mutate(Timepoint = replace(Timepoint, group_per_sample == "control", "HC")) %>%
        ungroup()

meta_sub <- meta_sub %>% 
        #mutate(age_numeric = as.numeric(as.factor(age))) %>%
        mutate(age_numeric = replace(age, age == "n/a", NA)) %>%
        mutate(age_numeric = as.numeric(sapply(strsplit(age_numeric, "_"), `[[`, 1)) + 2)

timepoint_dat <- meta_sub %>% select(sampleID, Timepoint, age_numeric)

meta <- left_join(meta, timepoint_dat)

meta <- meta %>% mutate(onset_group2 = NA) %>%
        mutate(onset_group2 = replace(onset_group2, days_after_onset < 7, "before7")) %>%
        mutate(onset_group2 = replace(onset_group2, days_after_onset == 0, "HC")) %>%
        mutate(onset_group2 = replace(onset_group2, days_after_onset >= 7 & days_after_onset < 15, "prejuncture7_15")) %>%
        mutate(onset_group2 = replace(onset_group2, days_after_onset >= 15 & days_after_onset <= 20, "juncture15_20")) %>%
        mutate(onset_group2 = replace(onset_group2, days_after_onset > 20, "after20"))

meta <- meta %>% mutate(severity_onsetgroup2 = paste(group_per_sample, onset_group2, sep = "_"))

meta %>% select(c("severity_onsetgroup2", keep_cols)) %>% distinct() %>%
        pull(severity_onsetgroup2) %>% table()

obj$Timepoint <- meta$Timepoint
obj$age_numeric <- meta$age_numeric
#obj$onset_group <- meta$onset_group
obj$onset_group2 <- meta$onset_group2
obj$severity_onsetgroup2 <- meta$severity_onsetgroup2
obj$who_per_sample <- as.numeric(obj$who_per_sample)

obj$group_per_sample <- relevel(factor(obj$group_per_sample), ref = "mild")

obj$id.celltype <- gsub(": ", "_", obj$id.celltype)

saveRDS(obj, OUT_PATH, compress = FALSE)
