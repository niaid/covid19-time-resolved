library(edgeR)

DGELIST_IN_PATH <- snakemake@input[[1]]
DGELIST_OUT_PATH <- snakemake@output[[1]]
pseudobulk_list <- readRDS(DGELIST_IN_PATH)


LIBSIZE_FILTER <- as.numeric(snakemake@params[["libsize_filter"]])
MIN_CELLS_PER_POOL <- snakemake@params[["min_cells_per_pool"]]
MIN_SAMPLES_PER_CELLTYPE <- snakemake@params[["min_samples_per_celltype"]] 


NormalizePseudobulk = function(dgelist, normalization_method = "RLE", 
                               group_col = NULL, minimum_gene_count = 1) {

  if(is.null(group_col)){
    group <- NULL
  }else{
    group <- dgelist$samples[[group_col]]
  }

  keep_genes = filterByExpr(dgelist, 
                            group=group,
                            min.count = minimum_gene_count,
                            min.total.count = 0)
  print(paste("keeping", sum(keep_genes), "genes"))
  dgelist = dgelist[keep_genes, keep.lib.sizes=FALSE]
  dgelist = calcNormFactors(dgelist, method = "RLE")
}



pseudobulk_list <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > LIBSIZE_FILTER & dge$samples$n_barcode > MIN_CELLS_PER_POOL]
})

n_samples_per_celltype <- sapply(pseudobulk_list, ncol)

pseudobulk_list <- pseudobulk_list[n_samples_per_celltype > MIN_SAMPLES_PER_CELLTYPE]

pseudobulk_list <- lapply(pseudobulk_list, NormalizePseudobulk)

saveRDS(pseudobulk_list, DGELIST_OUT_PATH)
