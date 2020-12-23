library(dplyr)
library(purrr)
library(edgeR)

#This is the primary function that you should use to do pooling

getPseudobulkList <- function(mat, celltypes, samples, pooling_function = c("sum", "mean"), 
                              meta = NULL, barcode_col_name = NULL, sample_level_meta_cols = NULL,
                              cell_level_meta_cols = NULL,
                              min_cells_per_pool = 10, 
                              min_samples_per_celltype = 10,
                              output_type = c("DGEList", "list")){

  cell_barcode_list <- getCellSampleBarcodesList(barcodes = colnames(mat), celltypes = celltypes, samples = samples,
                                           min_cells_per_pool = min_cells_per_pool, 
                                           min_samples_per_celltype = min_samples_per_celltype)
  
  cell_sample_nested_list <- poolCellSampleBarcodeList(mat, cell_barcode_list,
                                                       pooling_function = pooling_function,
                                                       meta = meta, barcode_col_name = barcode_col_name,
                                                       sample_level_meta_cols = sample_level_meta_cols,
                                                       cell_level_meta_cols = cell_level_meta_cols,
                                                       output_type = output_type)
  return(cell_sample_nested_list)

}

# The rest of these functions are called internally by getPseudobulkList, but can be used on their own


#Takes nested list barcodes (Celltype X Sample) from getCellSampleBarcodesList and returns nested pseudobul DGElists (Celltype X sample)
poolCellSampleBarcodeList <- function(mat, cell_sample_barcode_list, pooling_function = c("sum", "mean"), 
                                      meta = NULL, barcode_col_name = NULL, sample_level_meta_cols = NULL,
                                      cell_level_meta_cols = NULL,
                                      output_type = c("DGEList", "list")){

  cell_sample_nested_list <- lapply(cell_sample_barcode_list, 
                                    poolSampleBarcodeList, 
                                    mat = mat, pooling_function = pooling_function,
                                    meta = meta, barcode_col_name = barcode_col_name, 
                                    sample_level_meta_cols = sample_level_meta_cols, 
                                    cell_level_meta_cols = cell_level_meta_cols,
                                    output_type = output_type)
  return(cell_sample_nested_list)
}

# returns a nested list of cell barcodes for a given pool
# can filter such that only pools with a suffificent number of cells is maintained

getCellSampleBarcodesList <- function(barcodes, celltypes, samples, min_cells_per_pool = 10, min_samples_per_celltype = 10){
  
  cell_count_dat <- as.data.frame(table(celltypes = celltypes, samples = samples)) %>%
          mutate(n_cells_in_pool = Freq)
  cell_count_dat <- cell_count_dat %>% filter(n_cells_in_pool > min_cells_per_pool)

  sample_count_dat <- cell_count_dat %>%
          group_by(celltypes) %>%
          summarise(n_samples = n())

  keep_cell_types <- sample_count_dat %>% filter(n_samples > min_samples_per_celltype) %>%
          pull(celltypes)

  out_list <- list()
  for(celltype in keep_cell_types){
    out_list[[celltype]] <- list()
    keep_samples <- cell_count_dat %>% filter(celltypes == celltype) %>% pull(samples)
    for(SAMPLE in keep_samples){
      out_list[[celltype]][[SAMPLE]] <- barcodes[celltypes == celltype & samples == SAMPLE]
    }
  }

  return(out_list)
}
#Given a names list of cell barcodes, returns a gene X sample pseudobulk expression matrix
poolExpressionSampleBarcodeList <- function(mat, sample_barcode_list, pooling_function = c("sum", "mean")){
  pool_mat <- lapply(sample_barcode_list, function(barcodes){
    if(pooling_function == "sum"){
      Matrix::rowSums(mat[, barcodes])
    }else if(pooling_function == "mean"){
      Matrix::rowMeans(mat[, barcodes])
    }else{
      stop("Invalid pooling function: must be 'sum' or 'mean'")
    }
  })
  pool_mat <- do.call(cbind, pool_mat)
  return(pool_mat)
}

#Given a names list of cell barcodes, returns a dataframe of sample level metadata 
poolSampleMetadataSampleBarcodeList <- function(meta, sample_barcode_list, barcode_col_name, sample_level_meta_cols, keep_barcodes_sample_meta = TRUE){
  sample_names <- names(sample_barcode_list)
  meta_sub <- lapply(sample_names, function(sample_name){
          barcodes <- sample_barcode_list[[sample_name]]
          keep_rows <- meta[[barcode_col_name]] %in% barcodes
          meta_single <- meta[keep_rows, sample_level_meta_cols, drop = FALSE]
          meta_single <- distinct(meta_single)

          if(nrow(meta_single) > 1){
            stop(paste("The metadata in sample_level_meta_cols has multiple values for sample", sample_name))
          }
          if(keep_barcodes_sample_meta){
            meta_single[["barcodes"]] <- paste(barcodes, collapse = ",")
          }
          return(meta_single)
  })
  meta_sub <- do.call(rbind, meta_sub)
  rownames(meta_sub) <- sample_names
  return(meta_sub)
}

#Given a names list of cell barcodes, returns a nested list (Feature X Sample) of cell level metadata into a nested list
poolCellMetadataSampleBarcodeList <- function(meta, sample_barcode_list, barcode_col_name, cell_level_meta_cols){
  cell_meta <- lapply(sample_barcode_list, function(barcodes){
          keep_rows <- meta[[barcode_col_name]] %in% barcodes
          meta_sub <- meta[keep_rows, cell_level_meta_cols, drop = FALSE]
          as.list(meta_sub)
  })

  #make is so that 
  cell_meta <- purrr::transpose(cell_meta)
  return(cell_meta)
}


#Given a names list of cell barcodes, returns a list or DGEList with expression, sample level metadata and cell level metadata
poolSampleBarcodeList <- function(mat, sample_barcode_list, pooling_function = c("sum", "mean"), 
                                  meta = NULL, barcode_col_name = NULL, 
                                  sample_level_meta_cols = NULL, 
                                  cell_level_meta_cols = NULL, 
                                  keep_barcodes_sample_meta = TRUE,
                                  output_type = c("DGEList", "list")){

  pool_mat <- poolExpressionSampleBarcodeList(mat = mat, sample_barcode_list = sample_barcode_list,
                                                    pooling_function = pooling_function)
  
  n_barcodes_per_pool <- sapply(sample_barcode_list, length)
  if(!is.null(meta) & !is.null(sample_level_meta_cols)){
    sample_meta <- poolSampleMetadataSampleBarcodeList(meta = meta, sample_barcode_list = sample_barcode_list, 
                                                           barcode_col_name = barcode_col_name, 
                                                           sample_level_meta_cols = sample_level_meta_cols,
                                                           keep_barcodes_sample_meta = keep_barcodes_sample_meta)
    sample_meta$n_barcodes <- n_barcodes_per_pool
  
  }else{
    sample_meta <- data.frame(n_barcdodes = n_barcodes_per_pool)
  }

  if(!is.null(meta) & !is.null(cell_level_meta_cols)){
    cell_meta <- poolCellMetadataSampleBarcodeList(meta = meta, sample_barcode_list = sample_barcode_list, 
                                                   barcode_col_name = barcode_col_name, 
                                                   cell_level_meta_cols = cell_level_meta_cols)
  }else{
    cell_meta <- NULL
  }
  

  if(output_type == "DGEList"){
    out_list <- DGEList(counts = pool_mat, samples = sample_meta)
    #This is kind of a hack but it works for now
    out_list$cell_meta <- cell_meta
  }else if(output_type == "list"){
    out_list <- list(counts = pool_mat, samples = sample_meta, cell_meta = cell_meta)
  }else if(output_type == "ExpressionSet"){
    require(Biobase)
    out_list <- ExpressionSet(assayData = pool_mat, 
                              phenoData = AnnotatedDataFrame(data = sample_meta))
  }else{
    stop("Invalid output_type selected")
  }

  return(out_list)

}

