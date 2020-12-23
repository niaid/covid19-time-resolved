library(tidyverse)
library(edgeR)
library(Biobase)
library(reshape2)

### GSVA scores from filtered samples ########################################################
# read in covid all timepoints + healthy samples metadata
# input data from pseudobulk eset object -- output from pseudobulk.normalize.esetlist.R script
eset_list_all_WCTcourse <- readRDS("output/dge_lists/pbulk_eset_list_normalized_WCTcourse_metafiltered.rds")
eset_list_all_specificgating <- readRDS("output/dge_lists/pbulk_eset_list_normalized_specific_gating_metafiltered.rds")
names(eset_list_all_specificgating) <- paste("gated", names(eset_list_all_specificgating), sep = " ")
eset_list_all <- c(eset_list_all_WCTcourse, eset_list_all_specificgating)


### utility function #########################################################################
df_to_gsva_elist <- function(gsva_df, meta_eset){
  gsva_list <- group_split(gsva_df, celltype, .keep = TRUE)
  names(gsva_list) <- sapply(gsva_list, function(x){unique(x$celltype)})
  
  gsva_list <- lapply(gsva_list, function(dge){
    mat <- dcast(dge, pathway ~ sample, value.var = "module.score")
    rownames(mat) <- mat$pathway
    mat2 <- as.matrix(mat[,-1])
    meta <- meta_eset[[as.character(unique(dge$celltype))]]@phenoData@data
    meta <- meta[intersect(rownames(meta), unique(dge$sample)),]
    meta$Subject <- sapply(str_split(rownames(meta), "_"), function(x)x[1])
    mat2 <- mat2[,rownames(meta)]
    ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(meta))
  })
  
  # tidy meta data
  gsva_list <- lapply(gsva_list, function(dge){
    dge@phenoData@data$PC1_cat <- as.character(dge@phenoData@data$PC1_cat)
    dge@phenoData@data$PC1_cat <- replace(dge@phenoData@data$PC1_cat, 
                                          dge@phenoData@data$Class == "HC", "HC")
    dge@phenoData@data$PC1_cat <- factor(dge@phenoData@data$PC1_cat, levels = c("HC","PC1_low","PC1_high"))
    filter <- colnames(dge)
    dge <- dge
  })
  
  
  return(gsva_list)
}





### GSVA scores using filtered samples #####################################################
# the input is the gsva score tables
PC1_gsva <- readRDS("input/PC1_PC1_module_score_gsva_filtered_samples_genes.rds")
PC1mid_gsva <- readRDS("input/PC1_onset_group_interaction_PC1High-low_in_mid_module_score_gsva_filtered_samples_genes.rds")
PC1_onset_union_gsva <- readRDS("input/union_PC1_days_since_onset_module_score_gsva_filtered_samples_genes.rds")

PC1_gsva_esetlist <- df_to_gsva_elist(gsva_df = PC1_gsva, meta_eset = eset_list_all)
PC1mid_gsva_esetlist <- df_to_gsva_elist(gsva_df = PC1mid_gsva, meta_eset = eset_list_all)
PC1_onset_union_gsva_esetlist <- df_to_gsva_elist(gsva_df = PC1_onset_union_gsva, meta_eset = eset_list_all)

dir.create("output/module_score_gsva")
saveRDS(PC1_gsva_esetlist, "output/module_score_gsva/PC1_PC1_module_score_gsva_filtered_samples_genes_esetlist.rds")
saveRDS(PC1mid_gsva_esetlist, "output/module_score_gsva/PC1High-low_in_mid_module_score_gsva_filtered_samples_genes_esetlist.rds")
saveRDS(PC1_onset_union_gsva_esetlist, "output/module_score_gsva/union_PC1_days_since_onset_module_score_gsva_filtered_samples_genes_esetlist.rds")


