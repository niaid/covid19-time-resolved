library(edgeR)
library(Seurat)
library(matrixStats)
library(Matrix)
library(tidyverse)

#source
source("pseudobulk_pooling_functions.R")

SEURAT_IN_PATH <- "~/sg_data/users/martinsaj/COVID19/CITE5P_R_analysis/B1merge.Clustered.20200507.Rds"


MIN_CELLS_PER_POOL <- 10

#read in data
seurat_obj <- readRDS(SEURAT_IN_PATH)

meta <- seurat_obj@meta.data
rna <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
rm(seurat_obj)

samples <- paste(meta$Donor, meta$Timepoint, sep = "_")

#testing main function
pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta$mergedcelltype, 
                                     meta = meta,
                                     samples = samples, 
                                     min_cells_per_pool = 10, 
                                     min_samples_per_celltype = 10,
                                     barcode_col_name = "NewBarcode", 
                                     sample_level_meta_cols = c("Donor", "Age"),
                                     cell_level_meta_cols = NULL)

pseudobulk_list_filtered_libsize <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > 5e4]
  
})

sapply(pseudobulk_list_filtered_libsize, ncol)

library(txtplot)

for(nm in names(pseudobulk_list)){
  print(nm)
  txtdensity(log10(pseudobulk_list[[nm]]$samples$lib.size))
}

#testing with cell meta cols
pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta$mergedcelltype, 
                                     meta = meta,
                                     samples = samples, 
                                     min_cells_per_pool = 10, 
                                     min_samples_per_celltype = 10,
                                     barcode_col_name = "NewBarcode", 
                                     sample_level_meta_cols = c("Donor", "Age"),
                                     cell_level_meta_cols = c("Donor", "Age"))



#Testing component function

tmp <- getCellSampleBarcodesList(barcodes = meta$NewBarcode, celltypes = meta$mergedcelltype, samples = samples)

e <- poolExpressionSampleBarcodeList(rna, sample_barcode_list = tmp[[1]], pooling_function = "sum")

sample_meta <- 
        poolSampleMetadataSampleBarcodeList(meta = meta, sample_barcode_list = tmp[[1]], 
                                            barcode_col_name = "NewBarcode", 
                                            sample_level_meta_cols = c("Donor", "Age"))
cell_meta <- 
        poolCellMetadataSampleBarcodeList(meta = meta, sample_barcode_list = tmp[[1]], 
                                            barcode_col_name = "NewBarcode", 
                                            cell_level_meta_cols = c("Donor", "Age"))

test <- poolSampleBarcodeList(mat = rna, sample_barcode_list = tmp[[1]], 
                                          pooling_function = "sum", 
                                          meta = meta, barcode_col_name = "NewBarcode", 
                                          sample_level_meta_cols = "Age", 
                                          cell_level_meta_cols = "Age", 
                                          output_type = "DGEList")

test1 <- poolSampleBarcodeList(mat = rna, sample_barcode_list = tmp[[1]], 
                                          pooling_function = "sum")


test2 <- poolCellSampleBarcodeList(mat = rna, cell_sample_barcode_list = tmp, 
                                   pooling_function = "sum", 
                                   meta = meta, barcode_col_name = "NewBarcode", 
                                   sample_level_meta_cols = "Age",
                                   cell_level_meta_cols = "Age",
                                   output_type = "DGEList")

sample_meta_w_cell_meta <- sample_meta
sample_meta$asdf <- cell_meta$Donor

dge <- DGEList(counts =  - e, samples = sample_meta)
dge$cell_meta <- cell_meta

dge[1, ]
dge[, 1]
