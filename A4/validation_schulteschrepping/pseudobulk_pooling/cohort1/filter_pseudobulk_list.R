library(edgeR)

LIBSIZE_FILTER <- as.numeric(snakemake@params[["libsize_filter"]])
MIN_CELLS_PER_POOL <- snakemake@params[["min_cells_per_pool"]]
MIN_SAMPLES_PER_CELLTYPE <- snakemake@params[["min_samples_per_celltype"]] 

DGELIST_IN_PATH <- snakemake@input[[1]]

DGELISTS_OUT_DIR <- snakemake@output[[1]]

pseudobulk_list <- readRDS(DGELIST_IN_PATH)

pseudobulk_list <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > LIBSIZE_FILTER & dge$samples$n_barcode > MIN_CELLS_PER_POOL]
})

n_samples_per_celltype <- sapply(pseudobulk_list, ncol)

pseudobulk_list <- pseudobulk_list[n_samples_per_celltype > MIN_SAMPLES_PER_CELLTYPE]

for(nm in names(pseudobulk_list)){
  out_path <- paste0(DGELISTS_OUT_DIR, nm, "-pseudobulk_dgelist.rds")
  saveRDS(pseudobulk_list[[nm]], out_path)
}
