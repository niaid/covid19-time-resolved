library(tidyverse)
library(Biobase)


DAT_IN_PATH <- snakemake@input[[1]]
META_IN_PATH <- snakemake@input[["meta"]]
ESET_OUT_PATH <- snakemake@output[[1]]

meta <- readRDS(META_IN_PATH)
dat <- read_csv(DAT_IN_PATH)

dat <- dat %>% rename(batch_subj_time = Var1) %>%
        select(-c(X1, Timepoint, sample_id, Batch, days_since_symptoms_onset, severity))


mat <-  dat %>% 
        `rownames<-`(.$batch_subj_time) %>%
        select(-batch_subj_time) %>%
        as.matrix() %>%
        t()

colnames(mat) <- gsub("_PBMC", "", colnames(mat))

sample_cols <- c("Donor", "Age", "condition", "sex", 
               "ever_admitted_to_icu", "intubation_vent",
               "outcome", "date_of_onset_of_symptoms", "hospital_admission_date", 
               "hospital_discharge_date", "day_from_symptom_onset_to_hospitalization", 
               "number_of_days_hospitalized", "days_since_onset", "days_since_hospitalized", 
               "T0_crp", "batch", "severity", "outcome", "intubation_vent")

sample_cols <- c(sample_cols, "severity.outcome", "PC1", "PC2", "PLS1", 
                 "ANC.ALC.Ratio_.ratio", "D.Dimer_ng.mL", "IL.6",
                 "IP.10","Lactate.Dehydrogenase..LDH._U.L",
                 "Lymphocytes_x10.3.uL", "TNF.b",
                 "SpO2.FiO2_.ratio", "initial.Nucleocapsid",
                 "initial.Spike", "initial.days_from_time_point",
                 "end.Nucleocapsid", "end.Spike",
                 "end.days_from_discharge_or_death",
                 "spike.dicharge.delta.t0" ,
                 "nucleocapsid.dicharge.delta.t0" )

sample_cols <- c(sample_cols, "PC1_cat", "PC2_cat", "PLS1_cat")
sample_cols <- c(sample_cols, "Class", "cond_group", "Gender", "date_drawn", "Timepoint", "Batch")

meta_sub <- meta %>%
        select(sample_cols) %>%
        distinct()
meta_sub <- meta_sub %>% mutate(batch_subj_time = paste(Batch, Donor, Timepoint, sep = "_"))
rownames(meta_sub) <- meta_sub$batch_subj_time

stopifnot(sum(duplicated(meta_sub$batch_subj_time)) == 0)

pdat <- meta_sub[match(colnames(mat), meta_sub$batch_subj_time), ]
eset <- ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(pdat))

saveRDS(eset, ESET_OUT_PATH)
