library(tidyverse)
library(pheatmap)
library(edgeR)
library(circlize)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))

FIG_OUT_PATH <- "plots/CITE5p/all_batches/tcr/exhaustion/FIG4S.2020_08_21_exhaustion_signature_clustering_pbulk_updn_anno.pdf"

IN_PATHS <-
        c(CD8_Mem_sorted = "data/CITE5p/all_batches/differential_expression/2020_07_24/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Sorted-WCTcoursecelltype/CD8_Mem-pseudobulk_dgelist_normalized.rds",
         CD8_Mem_unsorted = "data/CITE5p/all_batches/differential_expression/2020_07_24/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Unsorted-WCTcoursecelltype/CD8_Mem-pseudobulk_dgelist_normalized.rds",
        CD8_Mem_sorted_expanded = "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Sorted-WCTcoursecelltype/CD8_Mem_expandedTRUE-pseudobulk_dgelist_normalized.rds",
        CD8_Mem_sorted_unexpanded = "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Sorted-WCTcoursecelltype/CD8_Mem_expandedFALSE-pseudobulk_dgelist_normalized.rds",
        CD8_Mem_unsorted_expanded = "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Unsorted-WCTcoursecelltype/CD8_Mem_expandedTRUE-pseudobulk_dgelist_normalized.rds",
        CD8_Mem_unsorted_unexpanded = "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups/t0_plus_healthy/pseudobulk_dgelists_normalized/Unsorted-WCTcoursecelltype/CD8_Mem_expandedFALSE-pseudobulk_dgelist_normalized.rds")

genesets <- readRDS("genesets/kegg_go_btm_reactome_foointerferon.rds")

wherry_exhaustion_dn_dat <- data.frame(genes = genesets$wherry_2007_exhaustion_dn, direction = "down", stringsAsFactors = F)
wherry_exhaustion_up_dat <- data.frame(genes = genesets$wherry_2007_exhaustion_up, direction = "up", stringsAsFactors = F)
wherry_exhaustion_dat <- bind_rows(wherry_exhaustion_dn_dat, wherry_exhaustion_up_dat) %>% distinct()

mckinney_cd2_exhaustion_up_dat <- data.frame(genes= genesets$mckinney_2015_cd8_cd2dn_exhaustion_up, direction = "up", stringsAsFactors = F)
mckinney_cd2_exhaustion_dn_dat <- data.frame(genes= genesets$mckinney_2015_cd8_cd2up_exhaustion_dn, direction = "down", stringsAsFactors = F)
mckinney_exhaustion_dat <- bind_rows(mckinney_cd2_exhaustion_dn_dat, mckinney_cd2_exhaustion_up_dat) %>% distinct()


#exhaustion_sigs$wherry_exhaustion_up_dn_rm_RPS_RPL <- 
#        exhaustion_sigs$wherry_exhaustion_up_dn[!grepl("RPS", exhaustion_sigs$wherry_exhaustion_up_dn) & ! grepl("RPL", exhaustion_sigs$wherry_exhaustion_up_dn)]

pbulk_list <- lapply(IN_PATHS, readRDS)

cpm_list <- lapply(pbulk_list, function(dge){
  dge <- dge[, !grepl("CHI", dge$sample$Donor)]
  cpm_mat <- cpm(dge, log = TRUE)                                
  list(cpm_mat = cpm_mat, samples = dge$samples)
})

exhaustion_sigs <- list(wherry_exhaustion_up_dn = wherry_exhaustion_dat,
                        mckinney_cd2_exhaustion_up_dn = mckinney_exhaustion_dat)

pdf(FIG_OUT_PATH)
for(sig in "wherry_exhaustion_up_dn"){
  for(cell_group in "CD8_Mem_sorted_expanded"){
    obj <- cpm_list[[cell_group]]
    sig_dat <- exhaustion_sigs[[sig]]
    keep_genes <- sig_dat[["genes"]]
    cpm_mat <- obj$cpm_mat
    cpm_mat <- cpm_mat[rownames(cpm_mat) %in% keep_genes, ]

    cpm_mat <- t(scale(t(cpm_mat)))

    cpm_mat <- removeBatchEffect(cpm_mat, batch = factor(obj$samples$Batch))

    if(nrow(cpm_mat) < 50){
      SHOW_ROWNAMES <- TRUE
    }else{
      SHOW_ROWNAMES <- FALSE
    }
    
    sample_annotation <- data.frame(Class = obj$samples$Class, Batch = obj$samples$Batch,
                             PC1_class = obj$samples$PC1_cat,
                             Age = obj$samples$Age)
    rownames(sample_annotation) <- colnames(cpm_mat)
    gene_annotation <- data.frame(direction = sig_dat$direction)
    rownames(gene_annotation) <- sig_dat$gene

    pheatmap(cpm_mat, main = paste(sig, cell_group, sep = "\n"), 
             show_rownames = SHOW_ROWNAMES, 
             annotation_col = sample_annotation,
             annotation_row = gene_annotation,
             color = col_fun(seq(from = -3, to = 3, length.out = 100)))
  }
}
dev.off()
#file.remove(FIG_OUT_PATH)

source("/hpcdata/sg/sg_data/users/rachmaninoffn2/dev/R/textboxplot_formula.R")

txtboxplot.formula(mpg~factor(cyl), data = mtcars)

for(cell_group in names(cpm_list)){
  obj <- cpm_list[[cell_group]]
  meta <- obj$samples
  print(cell_group)
  txtboxplot.formula(n_barcodes ~ Class, data = meta)
}

for(cell_group in names(cpm_list)){
  obj <- cpm_list[[cell_group]]
  meta <- obj$samples
  print(cell_group)
  txtboxplot.formula(lib.size ~ Class, data = meta)
}
