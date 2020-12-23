### utility funs ############################################################################
### get gsva score dataframe with meta from eset list
# For Schulte-Schrepping et al data

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
    if(nrow(dat) < 7){
      next()
    }
    
    dat$days_after_onset <- eset$days_after_onset
    dat$group <- eset$group
    dat$group_per_sample <- eset$group_per_sample
    dat$who_per_sample <- eset$who_per_sample
    dat$Timepoint <- eset$Timepoint
    dat$n_barcodes <- eset$n_barcodes
    dat$PC1 <- eset$PC1
    dat$age <- eset$age_numeric
    dat$donor <- eset$donor
    dat$sampleID <- eset$sampleID
    
    score_list[[cell]] <- dat
  }
  
  score_list_df <- bind_rows(score_list, .id = "celltype")
  return(score_list_df)
}

