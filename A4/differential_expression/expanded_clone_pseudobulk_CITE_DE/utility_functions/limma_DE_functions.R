library(limma)
library(Biobase)

RunLimma = function(eset, design_matrix, add_random_effect = FALSE, 
                        block_this_variable = NULL,
                        do_contrast_fit, my_contrast_matrix){ 
  require(limma)

  mat <- exprs(eset)
  meta <- pData(eset)
  #check that there are actually duplicates
  if(add_random_effect){
    n_duplicates <- sum(duplicated(meta[[block_this_variable]]))
    print(paste("There are", n_duplicates, "duplicates samples."))
    if(n_duplicates < 3){
      print(paste("Only", n_duplicates, "duplicates. Skipping duplicate correlation step"))
      add_random_effect <- FALSE
    }
  }

  if (add_random_effect == TRUE) {
    # voom call 1 get logcpm observational weights to feed into dup cor
    
    # Get covariance structure of subject repeat measurements (adding the random effects for subj) 
    cor = duplicateCorrelation(object = mat, 
                                    design = design_matrix, 
                                    block = meta[[block_this_variable]])
    
    # Fit model and contrasts 
    fit = lmFit(mat, design = design_matrix, 
                     block = meta[[block_this_variable]],
                     correlation = cor$consensus.correlation)
  } else { 
    # Fit model without random effect 
    fit = lmFit(mat, design = design_matrix)
  } 
  if (do_contrast_fit == TRUE) {
    c_fit = contrasts.fit(fit = fit, contrasts = my_contrast_matrix)
    eb_c_fit = eBayes(c_fit)
  } else { 
    eb_c_fit = eBayes(fit)
  }

  return(eb_c_fit)

}
