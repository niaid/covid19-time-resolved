library(Seurat)
library(edgeR)
library(tidyverse)

print(str(snakemake))

snakemake@source("pseudobulk_pooling_functions.R")

SEURAT_IN_PATH <- snakemake@input[[1]]

EXPANSION_DAT_IN_PATH <- snakemake@input[[2]]

DGELISTS_OUT_PATH <- snakemake@output[[1]]

#Sample filtering parameters
#KEEP_TIMEPOINTS <- snakemake@params[["keep_timepoints"]]
KEEP_SORTED <- snakemake@params[["keep_sorted"]]
KEEP_BATCHES <- snakemake@params[["keep_batches"]]
KEEP_CELLTYPES <- snakemake@params[["keep_celltypes"]]

#Parameters for pseudobulk pooling
LIBSIZE_FILTER <- as.numeric(snakemake@params[["libsize_filter"]])
MIN_CELLS_PER_POOL <- snakemake@params[["min_cells_per_pool"]]
MIN_SAMPLES_PER_CELLTYPE <- snakemake@params[["min_samples_per_celltype"]] 

#Cell annotation column
CELL_ANNOTATION_COLUMN <- snakemake@params[["cell_anno_col"]]

print("keep batches")
print(KEEP_BATCHES)
#print("keep timepoints")
#print(KEEP_TIMEPOINTS)
print("keep sorted")
print(KEEP_SORTED)

#read in data
expanded_dat <- read_tsv(EXPANSION_DAT_IN_PATH)
seurat_obj <- readRDS(SEURAT_IN_PATH)

#filter to just covid samples
#keep_cells <- WhichCells(seurat_obj, expression = Batch %in% KEEP_BATCHES & Timepoint %in% KEEP_TIMEPOINTS & Sorted == KEEP_SORTED)
meta_tmp <- seurat_obj@meta.data

expanded_dat <- expanded_dat %>% select(barcodeBatch, expanded)

meta_tmp <- left_join(expanded_dat, meta_tmp)
meta_tmp <- meta_tmp %>% filter(Batch %in% KEEP_BATCHES, Sorted == KEEP_SORTED, eval(parse(text = CELL_ANNOTATION_COLUMN)) %in% KEEP_CELLTYPES)
keep_cells <- WhichCells(seurat_obj, expression = barcodeBatch %in% meta_tmp$barcodeBatch)

seurat_obj <- subset(seurat_obj, cells = keep_cells)

seurat_obj$expanded <- meta_tmp$expanded[match(seurat_obj$barcodeBatch, meta_tmp$barcodeBatch)]

meta <- seurat_obj@meta.data

print("46")
rna <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
rm(seurat_obj)

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
                 "end.days_from_discharge_or_death")

sample_cols <- c(sample_cols, "PC1_cat", "PC2_cat", "PLS1_cat")
sample_cols <- c(sample_cols, "Class", "cond_group", "Batch", "Gender")
sample_cols <- c(sample_cols, "Timepoint")
sample_cols <- c(sample_cols, "onset_group")
sample_cols <- c(sample_cols, "severity.outcome2")
sample_cols <- c(sample_cols, "Mono_Classical")
sample_cols <- c(sample_cols, "PC1_onset_group")
sample_cols <- c(sample_cols, CELL_ANNOTATION_COLUMN)

sample_cols <- c(sample_cols, 
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
                 "spike.peak.Spike.corrected"
                 )

print("58")
print("sample_cols that aren't in metadata")
sample_cols[!sample_cols %in% colnames(meta)]

#This is what is different from the rest!!!

#Have to include batch because CHI014 is present in all of the batches
print("66")
samples <- paste(meta$Donor, meta$Timepoint, meta$Batch, sep = "_")

desired_pools <- paste0(meta[[CELL_ANNOTATION_COLUMN]], "_expanded", meta$expanded)

pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = desired_pools, 
                                     meta = meta,
                                     samples = samples, 
                                     min_cells_per_pool = MIN_CELLS_PER_POOL, 
                                     min_samples_per_celltype = MIN_SAMPLES_PER_CELLTYPE,
                                     barcode_col_name = "barcodeBatch", 
                                     sample_level_meta_cols = sample_cols,
                                     #cell_level_meta_cols = cell_cols, 
                                     pooling_function = "sum",
                                     output_type = "DGEList")

print("80")
pseudobulk_list <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > LIBSIZE_FILTER]
  
})

n_samples_per_celltype <- sapply(pseudobulk_list, ncol)
print("n_samples_per_celltype before final filtering")
print(n_samples_per_celltype)

pseudobulk_list <- pseudobulk_list[n_samples_per_celltype > MIN_SAMPLES_PER_CELLTYPE]

saveRDS(pseudobulk_list, DGELISTS_OUT_PATH)
