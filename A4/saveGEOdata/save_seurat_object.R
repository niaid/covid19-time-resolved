library(Seurat)
library(tidyverse)

SEURAT_IN_PATH <- "data/CITE5p/all_batches/2020_08_09.rmBuffSHD8.allcelltypelabels.merge.SNG.wmeta.WithinBatchClustered.Rds"

SEURAT_OUT_PATH <- "forSubmission/brescia_paper1_seurat.rds"

META_IN_PATH <- "forSubmission/covid19.metadata.paper1.RData"

obj <- readRDS(SEURAT_IN_PATH)

load(META_IN_PATH)

obj$sex <- NULL

rm_cols <- c("sex", "PC1class", "T0_crp", "t0_crp_date", "days_diff_draw_test", 
             "days_diff_pbmcdraw_T0crp", "test_date", "PC2", "PLS1", "IL.13", "TNF.b",
             "IL.6", "D.Dimer_ng.mL", "IP.10", 
             "Lymphocytes_x10.3.uL",
             "Lactate.Dehydrogenase..LDH._U.L", "ANC.ALC.Ratio_.ratio", "SpO2.FiO2_.ratio",
             "PC2_cat", "PLS1_cat", "severity.outcome2", "Mono_Classical",
             "date_of_onset_of_symptoms", "hospital_admission_date",
             "hospital_discharge_date",
             "date_of_death",
             "date_drawn",
             "date_received",
             "initial.test_date",
             "end.test_date",
             "nucleocapsid.peak.test_date",
             "spike.peak.test_date",
             "initial.Nucleocapsid",
             "initial.Spike",
             "initial.test_date",
             "initial.days_from_time_point",
             "end.Nucleocapsid",
             "end.Spike",
             "end.test_date",
             "end.days_from_discharge_or_death",
             "nucleocapsid.peak.Nucleocapsid",
             "nucleocapsid.peak.Spike",
             "nucleocapsid.peak.test_date",
             "nucleocapsid.peak.days_between_symptom_onset_and_test",
             "spike.peak.Nucleocapsid",
             "spike.peak.Spike",
             "spike.peak.test_date",
             "spike.peak.days_between_symptom_onset_and_test",
             "initial.days_from_symptom_onset",
             "initial.Nucleocapsid.corrected",
             "initial.Spike.corrected",
             "nucleocapsid.peak.Nucleocapsid.corrected",
             "spike.peak.Spike.corrected",
             "days_from_admission",
             "days_from_symptom_onset"
             )

for(nm in rm_cols){
        if(!nm %in% colnames(obj@meta.data)){
               # print(paste("warning: could not remove", nm "\nNot in @meta.data"))
               print(paste("warning: could not remove", nm, "\nNot in @meta.data"))
        }else{
                obj[[nm]] <- NULL 
        }
}

saveRDS()

meta_sub <- obj@meta.data %>%
        select(Donor, Timepoint, days_since_onset,  days_from_symptom_onset, days_since_hospitalized, days_from_admission) %>%
        distinct()

#meta_sub <- obj@meta.data %>%
#        select(Donor, Timepoint, days_since_onset, days_since_hospitalized) %>%
#        distinct()
#
williams_meta <- covid19.samples %>%
        filter(material_type == "PBMC") %>%
        select(subject_id, visit, days_from_symptom_onset_to_sample_drawn, days_from_admission_to_sample_drawn) %>%
        rename(Donor = subject_id, Timepoint = visit)

#
joined <- left_join(meta_sub, williams_meta)

joined %>% 
        select(Donor, Timepoint, days_since_onset, days_from_symptom_onset, 
               days_from_symptom_onset_to_sample_drawn) %>% 
        write_csv(path = "/hpcdata/sg/sg_data/users/rachmaninoffn2/scratch/covid_days_since_onset.csv")
#joined$days_since_onset == joined$days_from_symptom_onset_to_sample_drawn
#joined$days_since_hospitalized == joined$days_from_admission_to_sample_drawn
joined %>% 
        select(Donor, Timepoint, days_since_hospitalized, days_from_admission, 
               days_from_admission_to_sample_drawn) %>%
        write_csv(path = "/hpcdata/sg/sg_data/users/rachmaninoffn2/scratch/covid_days_since_admission.csv")

saveRDS(obj, SEURAT_OUT_PATH)

#obj@meta.data %>%
#        select(Donor, Timepoint, days_since_onset, days_since_hospitalized,
#               days_from_symptom_onset, days_from_admission) %>%
#head()
