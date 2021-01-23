library(Seurat)
library(edgeR)
library(tidyverse)

print(str(snakemake))

snakemake@source("pseudobulk_pooling_functions.R")

SEURAT_IN_PATH <- snakemake@input[[1]]

DGELISTS_OUT_PATH <- snakemake@output[[1]]

#Sample filtering parameters
#KEEP_TIMEPOINTS <- snakemake@params[["keep_timepoints"]]
KEEP_SORTED <- snakemake@params[["keep_sorted"]]
KEEP_BATCHES <- snakemake@params[["keep_batches"]]

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
seurat_obj <- readRDS(SEURAT_IN_PATH)

#filter to just covid samples
#keep_cells <- WhichCells(seurat_obj, expression = Batch %in% KEEP_BATCHES & Timepoint %in% KEEP_TIMEPOINTS & Sorted == KEEP_SORTED)
keep_cells <- WhichCells(seurat_obj, expression = Batch %in% KEEP_BATCHES &  Sorted == KEEP_SORTED)

print("39")

seurat_obj <- subset(seurat_obj, cells = keep_cells)

meta <- seurat_obj@meta.data

print("46")
rna <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
rm(seurat_obj)

sample_cols <- c("Donor", "Age", "condition", "sex", 
               "ever_admitted_to_icu", "intubation_vent",
               "outcome", "date_of_onset_of_symptoms", "hospital_admission_date", 
               "hospital_discharge_date", "day_from_symptom_onset_to_hospitalization", 
               "number_of_days_hospitalized", "days_since_onset", "days_since_hospitalized", 
               "T0_crp", "batch", "severity", "outcome", "intubation_vent", 
                 "PC1", "PC1_cat", "Class", "cond_group", "Batch", "Gender",
               "Timepoint", "onset_group", "PC1_onset_group")

sample_cols <- c(sample_cols, CELL_ANNOTATION_COLUMN)


print("58")
print("sample_cols that aren't in metadata")
sample_cols[!sample_cols %in% colnames(meta)]

#This is what is different from the rest!!!

#Have to include batch because CHI014 is present in all of the batches
print("66")
samples <- paste(meta$Donor, meta$Timepoint, meta$Batch, sep = "_")

pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta[[CELL_ANNOTATION_COLUMN]], 
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
