# plot PC1 and COVIDvsHC common LE genes
# use:
# subject.hm(exprs_mtx = Mono_Classical_mtx, meta = meta, celltype = "Mono_Classical", module = "IFN")
subject.hm <- function(exprs_mtx, meta, celltype, module){
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
  severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
  PC1.color <- colorRamp2(c(-3, 3.6), c("white", "#e97171"))
  timeonset.color <- colorRamp2(c(0, 40), c("white", "#6a2c70"))
  Class.color <- c("HC" = "#F8766D", "COVID" = "#eebb4d")
  PC1class.color <- c("HC"="#00BA38", "PC1_low"="#619CFF", "PC1_high"="#F8766D")
  
  # set how many genes to label, due to space limit
  if(nrow(exprs_mtx) < 20){
    nlabel = 10
  }else if (nrow(exprs_mtx) >= 20 & nrow(exprs_mtx) < 60){
    nlabel = 15
  }else {nlabel = round(nrow(exprs_mtx)/3)}
  label <- c(1:nlabel)
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    PC1class = meta$PC1_cat,
    PC1 = meta$PC1, 
    Severity_outcome = meta$severity.outcome2,
    col = list(PC1 = PC1.color,
               PC1class = PC1class.color,
               Severity_outcome = severity.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"), 
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$PC1_cat,
               col = col_fun,
               column_gap = unit(2, "mm")
  ) +
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 12)))
  return(p)
}


# plot all PC1 LE genes, and mark common LEs in both PC1 association and COVIDvsHC
# add pathway/genesets bars
# use:
# subject.hm.PC1(exprs_mtx = Mono_Classical_mtx_PC1, meta = meta, celltype = "Mono_Classical", module = "IFN",
#                PC1LEset = Mono_Classical_PC1_IFN_LE, covidLEset = Mono_Classical_covid_IFN_LE)
subject.hm.PC1 <- function(exprs_mtx, meta=meta, celltype, module, PC1LEset, covidLEset){
  # set marked genesets
  PC1_COVID <- rownames(exprs_mtx) %in% intersect(PC1LEset, covidLEset)
  names(PC1_COVID) <- rownames(exprs_mtx)
  
  # set how many genes to label, due to space limit
  if(nrow(exprs_mtx) < 20){
    nlabel = 10
  }else if (nrow(exprs_mtx) >= 20 & nrow(exprs_mtx) < 60){
    nlabel = 18
  }else {nlabel = round(nrow(exprs_mtx)/3)}  
  label <- which(rownames(exprs_mtx) %in% PC1LEset[1:nlabel])
  label_genes <- rownames(exprs_mtx)[label]
  
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
  severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
  PC1.color <- colorRamp2(c(-3, 3.6), c("white", "#e97171"))
  Class.color <- c("HC" = "#F8766D", "COVID" = "#eebb4d")
  PC1class.color <- c("HC"="#00BA38", "PC1_low"="#619CFF", "PC1_high"="#F8766D")
  
  # set heatmap annotation
  ha = HeatmapAnnotation(
    PC1class = meta$PC1_cat,
    PC1 = meta$PC1, 
    Severity_outcome = meta$severity.outcome2,
    col = list(PC1 = PC1.color,
               PC1class = PC1class.color,
               Severity_outcome = severity.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"), 
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$PC1_cat,
               row_split = PC1_COVID,
               col = col_fun,
               column_gap = unit(2, "mm")
  )+
    Heatmap(PC1_COVID + 0, name = "PC1-COVID-co", col = c("0" = "white", "1" = "blue"), 
            show_heatmap_legend = FALSE, width = unit(5, "mm")) +
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 12)))
  return(p)
}



### add row split for some genesets ##########################################################
### e.g. antigen presentation and for MHCI vs MHCII vs others
subject.hm.wsplit <- function(exprs_mtx, meta, celltype, module, rowsplit){
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
  severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
  PC1.color <- colorRamp2(c(-3, 3.6), c("white", "#e97171"))
  timeonset.color <- colorRamp2(c(0, 40), c("white", "#6a2c70"))
  Class.color <- c("HC" = "#F8766D", "COVID" = "#eebb4d")
  PC1class.color <- c("HC"="#00BA38", "PC1_low"="#619CFF", "PC1_high"="#F8766D")
  
  if(nrow(exprs_mtx) < 20){
    nlabel = 10
  }else if (nrow(exprs_mtx) >= 20 & nrow(exprs_mtx) < 60){
    nlabel = 18
  }else {nlabel = round(nrow(exprs_mtx)/3)}
  label <- c(1:nlabel)
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    PC1class = meta$PC1_cat,
    PC1 = meta$PC1, 
    Severity_outcome = meta$severity.outcome2,
    # Timeonset = meta$days_since_onset,
    # Class = meta$Class,
    col = list(PC1 = PC1.color,
               PC1class = PC1class.color,
               Severity_outcome = severity.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"), 
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$PC1_cat,
               col = col_fun,
               row_split = rowsplit,
               row_names_gp = gpar(fontsize = 8),
               column_gap = unit(2, "mm")
  )
    # rowAnnotation(foo = anno_mark(at = label,
    #                               labels = label_genes,
    #                               labels_gp = gpar(fontsize = 12)))
  return(p)
}



