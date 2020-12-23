### utility fun ############################################################################
### get gsva score dataframe with meta from eset list
getgsvascore_list_df <- function(gsva_esetlist, combined_genesets = combined_genesets){
  celltypes <- names(gsva_esetlist)
  score_list <- list()
  for(cell in celltypes){
    if(grepl("dblt", cell, ignore.case = T) | cell == "gated_out" | cell == "Unknown"){
      next()
    }
    print(cell)
    eset <- gsva_esetlist[[cell]]
    dat <- t(exprs(eset))
    dat <- as.data.frame(dat) %>% select(intersect(combined_genesets, colnames(dat)))
    # if(nrow(dat) < 7){
    #   next()
    # }
    
    dat$days_since_onset <- eset$days_since_onset
    dat$pc1_group <- eset$PC1_cat
    dat$severity_outcome <- eset$severity.outcome2
    dat$Timepoint <- eset$Timepoint
    dat$n_barcodes <- eset$n_barcodes
    dat$PC1 <- eset$PC1
    dat$Age <- eset$Age
    dat$Batch <- eset$Batch
    dat$Subject <- eset$Subject
    dat$steroid.use <- eset$steroid.use2
    
    dat <- dat %>% 
      filter(!is.na(pc1_group))
    score_list[[cell]] <- dat
  }
  
  score_list_df <- bind_rows(score_list, .id = "celltype")
  return(score_list_df)
}


