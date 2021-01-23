genesetFET <- function(hits,
                         all.sets, 
                         background, 
                         multiple.correct.method = p.adjust.methods,
                         rm_genes_not_in_background = TRUE,
                         min_set_size = 5
                         ){
  
  ### all.sets ###
  # a named list of gene sets
  # each list element is character vector of gene names or ID's
  
  ### background ###
  #a character vector of all gene names/ID's to which gene sets from all.sets are matched
  
  ### multiple.correct.method ###
  # the multiple correction method used. See help(p.adjust.methods)

  if(rm_genes_not_in_background){
    all.sets <- lapply(all.sets, function(x){x[x %in% background]})
  }

  all.sets <- all.sets[sapply(all.sets, length) > min_set_size]
  
  #dat_list <- lapply(all.sets, function(x) {singleEnrich(x, hits, background, rm_genes_not_in_background = FALSE)})

  dat_list <- vector(mode = "list", length = length(all.sets))
  names(dat_list) <- names(all.sets)
  #return(dat_list)

  for(gset_name in names(all.sets)){
    gset <- all.sets[[gset_name]]
    single_res <- singleEnrich(gset, hits, background, rm_genes_not_in_background = FALSE)
    if(is.null(single_res)){
      print("could not compute overlap for following geneset, no genes in background")
      print(gset_name)
    }
    dat_list[[gset_name]] <- single_res
    
  }

  dat_list <- dat_list[!is.null(dat_list)]

  combined_dat <- dplyr::bind_rows(dat_list, .id = "geneset_name")
  #combined_dat$geneset_name <- names(all.sets)

  combined_dat$p.adj <- p.adjust(combined_dat$P.value, method = multiple.correct.method)
  return(combined_dat)
}

singleEnrich <- function(gset, hits, background, return_hits = FALSE, return_overlap = TRUE, rm_genes_not_in_background = TRUE){
  if(any(!hits %in% background)){
    stop("Some hits not present in background")
  }
  # Performs a fisher's exact test for a single gene set
  # This is called by moduleFisher
  if(rm_genes_not_in_background){
    gset <- gset[gset %in% background]
  }
  in.mod <- background %in% hits
  
  in.gene.set <- background %in% gset

  if(sum(in.gene.set) == 0){
    print(paste("No genes in background, skipping geneset"))
    return(NULL)
  }
  
  contingency.tab <- table(in.gene.set, in.mod) #rows first argument in table, columns are second argument
  #print(contingency.tab)
  fisher_res <- fisher.test(contingency.tab, alternative = "greater")
  pval <- fisher_res$p.value
  odds_ratio <- fisher_res$estimate

  out_dat <- data.frame(P.value = pval, odds_ratio = odds_ratio, 
                        expected_overlap = contingency.tab[1,2] / sum(contingency.tab[1, ]) * length(hits),
                        observed_overlap = contingency.tab[2,2])
  if(return_hits){
    out_dat$hits <- paste(hits, collapse = ", ")
  }
  if(return_overlap){
    out_dat$overlap <- ifelse(out_dat$observed_overlap > 0, paste( background[in.mod & in.gene.set], collapse = ", "), NA)
  }
  return(out_dat)
}
