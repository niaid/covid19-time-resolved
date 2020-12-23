library(Seurat)
library(tidyverse)

SEURAT_IN_PATH <- "/hpcdata/sg/sg_data/users/liuc19/COVID19/CITE5p_Brescia_3batches/merge.SNG.wmeta.WithinBatchClustered.Rds"
META_WCRP_IN_PATH <- "data/metadata/cite5p_metadata_w_t0crp.tsv"
META_W_ENDPOINTS_IN_PATH <- "data/metadata/endpoints/2020_08_09/batch.sample.end.points.log10.RDS"


#WITHINCELLTYPE_CLUSTERS_IN_PATH <- "/hpcdata/sg/sg_data/users/martinsaj/COVID19/CITE5P_R_3Batches/withinBatchWithinCelltypeClustOut_200710/WCT.Labels.csv"

SEURAT_OUT_PATH <- "data/CITE5p/all_batches/2020_08_09.rmBuffSHD8.allcelltypelabels.merge.SNG.wmeta.WithinBatchClustered.Rds"

META_OUT_PATH <- "data/CITE5p/all_batches/2020_08_09.rmBuffSHD8.allcelltypelabels.merge.SNG.wmeta.WithinBatchClustered_metadata.Rds"

WCTCOURSE_PCNT_TOTAL_IN_PATH <- "/hpcdata/sg/sg_data/users/liuc19/COVID19/CITE5p_Brescia_3batches/auto_clusters/WCTcourse_ratio_to_all.20200730.csv"

source("scripts/nr/CITE5p/util/metadata_adding_functions.R")

seurat_obj <- readRDS(SEURAT_IN_PATH)

#remove buffy coat samples
keep_cells1 <- WhichCells(seurat_obj, expression = Subject %in% c("SHD8", "HDML_bc", "HDVO_bc"), invert = TRUE)
seurat_obj <- subset(seurat_obj, cells = keep_cells1)


wctcourse_dat <- read_csv(WCTCOURSE_PCNT_TOTAL_IN_PATH)
meta_w_crp <- read_tsv(META_WCRP_IN_PATH)
meta_w_endpoints <- readRDS(META_W_ENDPOINTS_IN_PATH)

meta_w_endpoints <- meta_w_endpoints %>%
        mutate(spike.peak.Spike = replace(spike.peak.Spike, is.na(spike.peak.Spike.corrected), NA)) %>%
        mutate(nucleocapsid.peak.Nucleocapsid= replace(nucleocapsid.peak.Nucleocapsid, 
                                                       is.na(nucleocapsid.peak.Nucleocapsid.corrected), NA))

#within_celltype_clusters_dat <- read_csv(WITHINCELLTYPE_CLUSTERS_IN_PATH)

meta_w_endpoints <- meta_w_endpoints %>%
        select(-c(date_of_onset_of_symptoms,hospital_admission_date, hospital_discharge_date, date_of_death))

#meta_w_endpoints <- meta_w_endpoints %>%
#        mutate(spike.discharge.delta.t0 = end.Spike - initial.Spike) %>%
#        mutate(nucleocapsid.discharge.delta.t0 = end.Nucleocapsid - initial.Nucleocapsid)

meta_w_endpoints <- meta_w_endpoints %>%
        mutate(PC1_cat = ifelse(PC1 > median(PC1), "PC1_high", "PC1_low")) %>%
        mutate(PC1_cat = relevel(factor(PC1_cat), "PC1_low")) %>%
        mutate(PC2_cat = ifelse(PC2 > median(PC2), "PC2_high", "PC2_low")) %>%
        mutate(PC2_cat = relevel(factor(PC2_cat), "PC2_low")) %>%
        mutate(PLS1_cat = ifelse(PLS1 > median(PLS1), "PLS1_high", "PLS1_low")) %>%
        mutate(PLS1_cat = relevel(factor(PLS1_cat), "PLS1_low"))


meta_w_crp <- meta_w_crp %>% rename(t0_crp_date = test_date)

keep_cols <- colnames(meta_w_endpoints)[!colnames(meta_w_endpoints) %in% colnames(meta_w_crp)]
meta_w_endpoints <- meta_w_endpoints[, c("subject_id", keep_cols)]

meta_w_crp <- left_join(meta_w_crp, meta_w_endpoints, by = "subject_id")

meta <- seurat_obj@meta.data

meta$barcodeBatch <- colnames(seurat_obj)

meta <- meta %>%
        mutate(WCTmergedcelltype = as.character(WCTmergedcelltype)) %>%
        mutate(WCTcoursecelltype= as.character(WCTcoursecelltype)) %>% 
        mutate(celltypeQC = as.logical(as.character(celltypeQC))) %>%
        mutate(WCTcoursecelltype = replace(WCTcoursecelltype, !celltypeQC, "gated_out")) %>%
        mutate(WCTmergedcelltype = replace(WCTmergedcelltype, !celltypeQC, "gated_out"))

meta <- joinMeta(cell_meta = meta,
                          sample_meta = meta_w_crp,
                          cell_meta_shared_colnames = c("Donor","Timepoint", "Age"),
                          sample_meta_shared_colnames = c("subject_id","visit",  "age"),
                          shared_colnames_output_names = c("Donor", "Timepoint", "Age"))

numeric_cols <- c("Age", "T0_crp", "days_since_onset", "days_since_hospitalized")
factor_cols <- c("condition", "Gender", "batch", "severity", "outcome", "intubation_vent")


for(nm in numeric_cols){
  if(!is.numeric(meta[[nm]])){
          print(nm)
    meta[[nm]] <- as.numeric(as.character(meta[[nm]]))
  }
}
for(nm in factor_cols){
        print(nm)
  meta[[nm]] <- droplevels(factor(meta[[nm]]))
}

# need to make new column that incorporates both hc and severity
meta$Class <- relevel(factor(meta$Class), "HC")
group <- as.character(meta$Class)
group[group == "COVID"] <- as.character(meta$severity)[group == "COVID"]
group <- relevel(factor(group), ref = "HC")
meta$cond_group <- group

#do releveling. 
meta$severity.outcome <- relevel(factor(meta$severity.outcome), "Severe-Alive")

#Add onset group
meta <- meta %>%
        mutate(onset_group = NA) %>%
        mutate(onset_group = replace(onset_group, days_since_onset < 17, "early")) %>%
        mutate(onset_group = replace(onset_group, days_since_onset >= 17 & days_since_onset <= 23, "mid")) %>%
        mutate(onset_group = replace(onset_group, days_since_onset > 23, "late"))

#meta %>% select(Donor, Timepoint, Batch, onset_group, days_since_onset) %>% distinct

#add severity.outcome column that has HC for healthy controls
meta <- meta %>%
        mutate(severity.outcome2 = as.character(severity.outcome)) %>%
        mutate(severity.outcome2 = replace(severity.outcome2, Class =="HC", "HC")) %>%
        mutate(severity.outcome2 = relevel(factor(severity.outcome2), "HC"))

#Add combined PC1 and onset group
meta <- meta %>%
        mutate(PC1_onset_group = paste(PC1_cat, onset_group, sep = "_"))

#meta %>% select(Donor, Timepoint, PC1_onset_group) %>%
#        distinct() %>%
#        pull(PC1_onset_group) %>%
#        table()

#add in WCTcourse monocytes pcnt total

wctcourse_dat <- wctcourse_dat %>% rename(Donor = Subject) %>%
        select(Donor, Timepoint, Batch, Mono_Classical)

meta_w_mono <- left_join(meta, wctcourse_dat)

head(colnames(seurat_obj) == meta_w_mono$barcodeBatch)
stopifnot(all(colnames(seurat_obj) == meta_w_mono$barcodeBatch))
for(nm in colnames(meta_w_mono)){
  seurat_obj[[nm]] <- meta_w_mono[[nm]]
}
head(colnames(seurat_obj) == seurat_obj$barcodeBatch)
stopifnot(all(colnames(seurat_obj) == seurat_obj$barcodeBatch))

saveRDS(meta_w_mono, META_OUT_PATH)
saveRDS(seurat_obj, SEURAT_OUT_PATH, compress = FALSE)
