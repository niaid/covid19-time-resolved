library(limma)
library(edgeR)

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

make_design_matrix <- function(dgelist, form){
  model.matrix(form, data = dgelist$samples)
}


RunVoomLimma = function(dgelist, design_matrix, add_random_effect = FALSE, 
                        block_this_variable = NULL,
                        do_contrast_fit, my_contrast_matrix){ 
  require(limma)
  require(edgeR)

  #check that there are actually duplicates
  if(add_random_effect){
    n_duplicates <- sum(duplicated(dgelist$samples[[block_this_variable]]))
    print(paste("There are", n_duplicates, "duplicates samples."))
    if(n_duplicates < 3){
      print(paste("Only", n_duplicates, "duplicates. Skipping duplicate correlation step"))
      add_random_effect <- FALSE
    }
  }

  if (add_random_effect == TRUE) { 
    # voom call 1 get logcpm observational weights to feed into dup cor
    v <- voom(counts = dgelist, normalize.method = "none",
                   design = design_matrix, save.plot = T, plot = F)
    
    # Get covariance structure of subject repeat measurements (adding the random effects for subj) 
    cor = duplicateCorrelation(object = v, 
                                    design = design_matrix, 
                                    block = dgelist$samples[[block_this_variable]])
    
    # rerun voom with the blocking variable and estimated correlation
    v = voom(dgelist, design = design_matrix, 
                  normalize.method = "none", 
                  block = dgelist$samples[[block_this_variable]], 
                  correlation=cor$consensus.correlation, plot = F, save.plot = T)
    
    # Fit model and contrasts 
    fit = lmFit(v, design = design_matrix, 
                     block = dgelist$samples[[block_this_variable]],
                     correlation = cor$consensus.correlation)
  } else { 
    v <- voom(counts = dgelist, normalize.method = "none",
                   design = design_matrix, save.plot = T, plot = F)
    # Fit model without random effect 
    fit = lmFit(v, design = design_matrix)
  } 
  if (do_contrast_fit == TRUE) {
    c_fit = contrasts.fit(fit = fit, contrasts = my_contrast_matrix)
    eb_c_fit = eBayes(c_fit)
  } else { 
    eb_c_fit = eBayes(fit)
  }

  return(eb_c_fit)

}
#  
# Usage 
# ebfit = RunVoomLimma(dgelist = d1sum_dge, design_matrix = d1m,
#   add_random_effect = T, block_this_variable = "subject", do_contrast_fit = T, 
#   my_contrast_matrix = contrast_matrix, dgelist$samples = d1)
# returns a list of linear model results. 



############################################################################################## 
##### Alternate pipeline with weights based on number of cells used to make pseudobulk library

# subset to celltypes used in model. 
# celltypes  = GetCelltypes(SeuratObject = h1, celltype = "celltype_joint")
# mytable = GetSubjectCelltypeTable(Seurat.Object = h1, celltype_column = "celltype_joint", sample_column = "sample")
# get sample weights as freq of cell total 
