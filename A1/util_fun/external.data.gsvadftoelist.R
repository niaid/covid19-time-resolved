### utility functions #########################################################################
### gsva dataframe to eset list with metadata
# For Schulte-Schrepping et al data
df_to_gsva_elist <- function(gsva_df, meta_eset){
  gsva_list <- group_split(gsva_df, celltype, .keep = TRUE)
  names(gsva_list) <- sapply(gsva_list, function(x){unique(x$celltype)})
  
  gsva_list <- lapply(gsva_list, function(dge){
    mat <- dcast(dge, pathway ~ sample, value.var = "module.score")
    rownames(mat) <- mat$pathway
    mat2 <- as.matrix(mat[,-1])
    meta <- meta_eset[[as.character(unique(dge$celltype))]]$samples
    meta <- meta[intersect(rownames(meta), unique(dge$sample)),]
    # meta$Subject <- sapply(str_split(rownames(meta), "_"), function(x)x[1])
    mat2 <- mat2[,rownames(meta)]
    stopifnot(identical(colnames(mat2),rownames(meta)))
    ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(meta))
  })
  
  return(gsva_list)
}





