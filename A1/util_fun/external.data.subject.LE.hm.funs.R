### pbulk heatmaps ##############################################################################
### utility functions ###########################################################################
library(ComplexHeatmap)
library(viridis)
library(circlize)

# For Schulte-Schrepping et al pbulk data
# plot all PC1 LE genes, and mark which are PC1 specific and which are also in COVIDvsHC
# add pathway/genesets bars
# use:
# subject.hm.PC1(exprs_mtx = Mono_Classical_mtx_PC1, meta = meta, celltype = "Mono_Classical", module = "IFN",
#                PC1LEset = Mono_Classical_PC1_IFN_LE, covidLEset = Mono_Classical_covid_IFN_LE)
subject.hm.PC1 <- function(exprs_mtx, meta=meta, celltype, module, 
                           PC1LEset, covidLEset, rowsplit=1, othergroup=NULL){
  # set marked genesets
  PC1highvslow <- rownames(exprs_mtx) %in% PC1LEset
  COVIDvsHC <- rownames(exprs_mtx) %in% covidLEset
  # othergroup <- rownames(exprs_mtx) %in% c(othergroup)
  names(COVIDvsHC) <- rownames(exprs_mtx)
  
  PC1_COVID <- rownames(exprs_mtx) %in% intersect(PC1LEset, covidLEset)
  names(PC1_COVID) <- rownames(exprs_mtx)
  
  if(rowsplit == 1){
    split = PC1_COVID
  } else {split = othergroup}  
  
  if(nrow(exprs_mtx) < 20){
    nlabel = 10
  }else if (nrow(exprs_mtx) >= 20 & nrow(exprs_mtx) < 60){
    nlabel = 18
  }else {nlabel = round(nrow(exprs_mtx)/3)}  
  label <- which(rownames(exprs_mtx) %in% union(PC1LEset, covidLEset)[1:nlabel])
  label_genes <- rownames(exprs_mtx)[label]
  
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
  who.color <- colorRamp2(c(0, 7), c("white", "#e97171"))
  severity.color <- c("control"="#00BA38", "mild"="#619CFF", "severe"="#F8766D")
  
  # set heatmap annotation
  ha = HeatmapAnnotation(
    severity = meta$group,
    who = meta$who_per_sample, 
    col = list(who = who.color,
               severity = severity.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"), 
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$group,
               row_split = split,
               col = col_fun,
               column_gap = unit(2, "mm"),
               row_gap = unit(2, "mm")
  )+
    Heatmap(PC1_COVID + 0, name = "severity-COVID-co", col = c("0" = "white", "1" = "blue"), 
            show_heatmap_legend = FALSE, width = unit(5, "mm")) +
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 12)))
  return(p)
}

