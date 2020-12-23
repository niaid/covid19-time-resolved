#this is tested using R 3.6.1 on a high-performance comupting node with 8 cores and at least 320 gb of ram. 
library("Seurat") #load Seurat 3.2.2
library("dplyr")
library("matrixStats")
library('tidyverse')
library('parallelDist')
library("pheatmap")
library("viridis")
library("ggridges")
library("magrittr")
library("genefilter")
library("DescTools")
library("scico")

#### Recluster on RNA within the T cell subsets ####

## import the seurat object downloaded from GEO
merge.SNG = readRDS("SeuratObjects/AllBatches_SeuratObj.rds")

## import fgsea results tables
t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable = read.table("Metadata/2020_08_09_t0_plus_healthy--healthy_vs_covid--COVID-Healthy.tsv", header = TRUE, sep = "\t")
t0_plus_t0_covid_only.PC1.fgseaatable = read.table("Metadata/2020_08_09_t0_covid_only--PC1--PC1.tsv", header = TRUE, sep = "\t")

### RNA clustering
library(genefilter)
f1 <- pOverA(0.01, 2)
ffun <- filterfun(f1)

dir.create("withinWCTCelltypeRNAclustPathLeadingEdge")
Pathway.leadingEdge = union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable$leadingEdge),  " "))),
                              unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$leadingEdge),  " ")))
                            )

#####  cluster RNA within the Total CD4

#filter out the cell subsets of interest
TotalCD4.LG = subset(merge.SNG, WCTcoursecelltype %in% c("CD4_Naive","CD4_Mem","Treg"))

#select just the leading edge genes for the pathway enrichment results of the relevant cell subsets
tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG = unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
  which(
    t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "btm_M4.1_cell cycle (I)",
      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
      "KEGG_RIBOSOME",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_OXIDATIVE_PHOSPHORYLATION",
      "reactome_Fatty acid metabolism"
    ) 
    &
      t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
        "CD4_Naive",
        "CD4_Mem",
        "Treg"
      )  
  ),]$leadingEdge),  " ")))
tO.PC1.SelectedPath.leadingEdge.LG = unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
  which(
    t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "btm_M4.1_cell cycle (I)",
      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
      "KEGG_RIBOSOME",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_OXIDATIVE_PHOSPHORYLATION",
      "reactome_Fatty acid metabolism"
    ) 
    &
      t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
        "CD4_Naive",
        "CD4_Mem",
        "Treg"
      )  
  ),]$leadingEdge),  " ")))

#run the clustering within the selected subset using filtered genes
TotalCD4.LG = ScaleData(TotalCD4.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG,tO.PC1.SelectedPath.leadingEdge.LG)), vars.to.regress=c("Donor","Batch"))
TotalCD4.LG <- RunPCA(TotalCD4.LG, assay="RNA", slot = "scale.data", 
                      features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG,tO.PC1.SelectedPath.leadingEdge.LG)), reduction.name = "pcaledge")
TotalCD4.LG <- FindNeighbors(TotalCD4.LG, dims = 1:15, reduction = "pcaledge")
TotalCD4.LG <- FindClusters(TotalCD4.LG, resolution = c(1), algorithm = 1)

#plottting average expression
source("utilityFunctions/AverageExpression_MeanOnly.r")
TotalCD4.LG.ADTfilt = genefilter(GetAssayData(TotalCD4.LG[["limmaCITE"]], slot = "data"), ffun)
TotalCD4.LGAver = AverageExpression_MeanOnly(TotalCD4.LG, return.seurat=F)
mat_breaks.temp <- quantile_breaks(as.matrix(TotalCD4.LGAver$limmaCITE[TotalCD4.LG.ADTfilt,]), n = 101)
ggsave(filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.WCTaverExp_withinCelltypeClust_RNA.1.pdf"), width = 8, height = 15,
       plot = pheatmap(TotalCD4.LGAver$limmaCITE[TotalCD4.LG.ADTfilt,], scale = "none", border_color=NA, viridis(n=100), breaks = mat_breaks.temp))
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$WCTmergedcelltype), margin = 2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPercentageperWCTmergedcelltypelabels.1.pdf"), width = 4, height = 6,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$adjustedcelltype), margin =  2)  %>% '*'(100) %>% round(2),
         filename =  paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeClustPercentageperAdjustedcelltypelabels.1.pdf"),width = 3, height = 6,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, paste(TotalCD4.LG$Donor, TotalCD4.LG$Timepoint, sep="_")), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPercentageperDonorTimepoint.1.pdf"), width = 10, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
TotalCD4.LG$WCTRNAclusters = paste0("TotalCD4.LG_",TotalCD4.LG$RNA_snn_res.1)
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$Class), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPerClass.1.pdf"), width = 2, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$severity.outcome), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPerSeverity.outcome.1.pdf"), width = 6, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$PC1_cat), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPerPC1_Cat.1.pdf"), width = 2, height = 12,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
write.csv(table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$sample_id),
          file = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.withinCelltypeRNAClustPerSample.csv"))

#Find the top gene markers of each cluster
TotalCD4.LG.RNAmarkers = FindAllMarkers(TotalCD4.LG, only.pos = TRUE)
write.csv(TotalCD4.LG.RNAmarkers, file = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.RNAMarkers_byRNAclust.csv"))
RNAcelltypeheatgenes = TotalCD4.LG.RNAmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
ggsave(filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD4.LG.COVID.batchaverExp_withinCelltypeClust_RNA.1.RNA.pdf"), width = 8, height = 16,
       plot = pheatmap(TotalCD4.LGAver$RNA[unique(RNAcelltypeheatgenes$gene),], scale = "row", border_color=NA, fontsize_row = 8))

#create Umap in the space of the selected genes
TotalCD4.LG <- RunUMAP(TotalCD4.LG, assay = "RNA", reduction = "pcaledge",
                    n.neighbors = 15, min.dist = 0.01, spread =  5, dims = 1:15)
TotalCD4.LGumap <- DimPlot(TotalCD4.LG, 
                        group.by = "RNA_snn_res.1", label = TRUE) + NoLegend()
TotalCD4.LGumap <- AugmentPlot(TotalCD4.LGumap, width = 4, height = 4, dpi = 300)
ggsave(TotalCD4.LGumap, filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.pdf", width = 4, height = 4)

# calculate module scores of the pathways in the cells and plot
TotalCD4.LGpaths <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      "HALLMARK_INFLAMMATORY_RESPONSE",
                      "btm_M4.1_cell cycle (I)",
                      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
                      "KEGG_RIBOSOME",
                      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                      "KEGG_OXIDATIVE_PHOSPHORYLATION",
                      "reactome_Fatty acid metabolism")
TotalCD4.LG <- AddModuleScore(TotalCD4.LG, features = genesets[TotalCD4.LGpaths], name = "pathway")
ggsave(DotPlot(TotalCD4.LG, features =  c("pathway1","pathway2","pathway3","pathway4","pathway5","pathway6","pathway7","pathway8")) + RotatedAxis(),
       filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.DotPlot.pdf", width = 4, height = 4
       )
ggsave(AugmentPlot(FeaturePlot(TotalCD4.LG, features = "pathway1", order = TRUE), width = 4, height = 4, dpi = 300), filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.HALLMARK_TNFA_SIGNALING_VIA_NFKB.pdf", width = 4, height = 4)
ggsave(AugmentPlot(FeaturePlot(TotalCD4.LG, features = "pathway2", order = TRUE), width = 4, height = 4, dpi = 300), filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.HALLMARK_INFLAMMATORY_RESPONSE.pdf", width = 4, height = 4)
ggsave(AugmentPlot(FeaturePlot(TotalCD4.LG, features = "pathway3", order = TRUE), width = 4, height = 4, dpi = 300), filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.btm_M4.1_cellcycleI.pdf", width = 4, height = 4)
ggsave(AugmentPlot(FeaturePlot(TotalCD4.LG, features = "pathway4", order = TRUE), width = 4, height = 4, dpi = 300), filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.GO_RESPONSE_TO_TYPE_I_INTERFERON.pdf", width = 4, height = 4)
ggsave(AugmentPlot(FeaturePlot(TotalCD4.LG, features = "pathway5", order = TRUE), width = 4, height = 4, dpi = 300), filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LGumap.WCTRNAclusters.KEGG_RIBOSOME.pdf", width = 4, height = 4)

FeaturePlot(TotalCD4.LG, features = "pathway2", order = TRUE)
FeaturePlot(TotalCD4.LG, features = "pathway3", order = TRUE)

ggsave(FeaturePlot(TotalCD4.LG, 
                   features = c("JUN"),
                   order = TRUE, ncol = 3),
       filename = "TotalCD4.LGheatplot.JUN.WCTRNAclusters.pdf", width = 6, height = 6)
ggsave(FeaturePlot(TotalCD4.LG, 
                   features = c("MX1"),
                   order = TRUE, ncol = 3),
       filename = "TotalCD4.LGheatplot.MX1.WCTRNAclusters.pdf", width = 6, height = 6)
ggsave(FeaturePlot(TotalCD4.LG, 
                   features = c("MKI67"),
                   order = TRUE, ncol = 3),
       filename = "TotalCD4.LGheatplot.MKI67.WCTRNAclusters.pdf", width = 6, height = 6)

TotalCD4.LG.UnSort = subset(TotalCD4.LG, Sorted == "N")
Idents(TotalCD4.LG.UnSort) <- "Class"
TotalCD4.LG.UnSort.30k <- SubsetData(TotalCD4.LG.UnSort, max.cells.per.ident = 30000)
TotalCD4.LGumapClass <- DimPlot(subset(TotalCD4.LG.UnSort.30k, Sorted == "N"), 
                             group.by = "RNA_snn_res.1", label = TRUE, split.by = "Class") + NoLegend()
ggsave(TotalCD4.LGumapClass, filename = "TotalCD4.LGumapClass.WCTRNAclusters.pdf", width = 12, height = 6)
TotalCD4.LGumapSeverityOutcome <- DimPlot(subset(TotalCD4.LG, Sorted == "N"), 
                                       group.by = "WCTRNAclusters", label = TRUE, split.by = "severity.outcome") + NoLegend()
ggsave(TotalCD4.LGumapSeverityOutcome, filename = "TotalCD4.LGumapSeverity.WCTRNAcelltype.pdf", width = 12, height = 6)
FeaturePlot(subset(TotalCD4.LG, Sorted == "N"), features = c("MKI67", "ISG15", "GZMA", "FOXP3"))
FeatureScatter(subset(TotalCD4.LG, Sorted =="N"), feature1 = "limmacite_CD25", "limmacite_CD127")
DotPlot(subset(TotalCD4.LG, Sorted == "N"), features = c("limmacite_CD25", "limmacite_CD127", "limmacite_HLA-DR",
                                                      "limmacite_CD71","limmacite_CD278","limmacite_CD69",
                                                      "limmacite_CD38", "limmacite_CD194", "limmacite_KLRG1",
                                                      "MKI67", "ISG15", "GZMA", "FOXP3", "TCF7", "TNFAIP3", "SELL", "CCR7"),
        group.by = "WCTRNAclusters") + RotatedAxis()


TotalCD4.LG.RNA.Labels = data.frame(WCTRNAclusters = TotalCD4.LG$WCTRNAclusters)
write.csv(TotalCD4.LG.RNA.Labels, file = "withinWCTCelltypeRNAclustPathLeadingEdge/RNAclustLabels/TotalCD4.LG.RNA.Labels.csv")


#####  cluster RNA within the Total CD8

#filter out the cell subsets of interest
TotalCD8.LG = subset(merge.SNG, WCTcoursecelltype %in% c("CD8_Naive","CD8_Mem"))

#select just the leading edge genes for the pathway enrichment results of the relevant cell subsets
tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG.CD8 = unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
  which(
    t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "btm_M4.1_cell cycle (I)",
      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
      "KEGG_RIBOSOME",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_OXIDATIVE_PHOSPHORYLATION",
      "reactome_Fatty acid metabolism"
    ) 
    &
      t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
        "CD8_Naive",
        "CD8_Mem"
      )  
  ),]$leadingEdge),  " ")))
tO.PC1.SelectedPath.leadingEdge.LG.CD8 = unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
  which(
    t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "btm_M4.1_cell cycle (I)",
      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
      "KEGG_RIBOSOME",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_OXIDATIVE_PHOSPHORYLATION",
      "reactome_Fatty acid metabolism"
    ) 
    &
      t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
        "CD8_Naive",
        "CD8_Mem"
      )  
  ),]$leadingEdge),  " ")))

#run the clustering within the selected subset using filtered genes
TotalCD8.LG = ScaleData(TotalCD8.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG.CD8,tO.PC1.SelectedPath.leadingEdge.LG.CD8)), vars.to.regress=c("Donor","Batch"))
TotalCD8.LG <- RunPCA(TotalCD8.LG, assay="RNA", slot = "scale.data",  features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG.CD8,tO.PC1.SelectedPath.leadingEdge.LG.CD8)), reduction.name = "pcaledge")
TotalCD8.LG <- FindNeighbors(TotalCD8.LG, dims = 1:15, reduction = "pcaledge")
TotalCD8.LG <- FindClusters(TotalCD8.LG, resolution = c(1), algorithm = 1)

#plottting average expression
TotalCD8.LG.ADTfilt = genefilter(GetAssayData(TotalCD8.LG[["limmaCITE"]], slot = "data"), ffun)
TotalCD8.LGAver = AverageExpression_MeanOnly(TotalCD8.LG, return.seurat=F)
mat_breaks.temp <- quantile_breaks(as.matrix(TotalCD8.LGAver$limmaCITE[TotalCD8.LG.ADTfilt,]), n = 101)
ggsave(filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.WCTaverExp_withinCelltypeClust_RNA.1.pdf"), width = 8, height = 15,
       plot = pheatmap(TotalCD8.LGAver$limmaCITE[TotalCD8.LG.ADTfilt,], scale = "none", border_color=NA, viridis(n=100), breaks = mat_breaks.temp))
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$WCTmergedcelltype), margin = 2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPercentageperWCTmergedcelltypelabels.1.pdf"), width = 4, height = 6,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$adjustedcelltype), margin =  2)  %>% '*'(100) %>% round(2),
         filename =  paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeClustPercentageperAdjustedcelltypelabels.1.pdf"),width = 3, height = 6,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, paste(TotalCD8.LG$Donor, TotalCD8.LG$Timepoint, sep="_")), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPercentageperDonorTimepoint.1.pdf"), width = 10, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
TotalCD8.LG$WCTRNAclusters = paste0("TotalCD8.LG_",TotalCD8.LG$RNA_snn_res.1)
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$Class), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPerClass.1.pdf"), width = 2, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$severity.outcome), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPerSeverity.outcome.1.pdf"), width = 6, height = 10,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
pheatmap(prop.table(x = table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$PC1_cat), margin =  2)  %>% '*'(100) %>% round(2),
         filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPerPC1_Cat.1.pdf"), width = 2, height = 12,
         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
write.csv(table(TotalCD8.LG$RNA_snn_res.1, TotalCD8.LG$sample_id),
    file = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.withinCelltypeRNAClustPerSample.csv"))

#Find the top gene markers of each cluster
TotalCD8.LG.RNAmarkers = FindAllMarkers(TotalCD8.LG, only.pos = TRUE)
write.csv(TotalCD8.LG.RNAmarkers, file = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.RNAMarkers_byRNAclust.csv"))
RNAcelltypeheatgenes = TotalCD8.LG.RNAmarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
ggsave(filename = paste0("withinWCTCelltypeRNAclustPathLeadingEdge/","TotalCD8.LG.COVID.batchaverExp_withinCelltypeClust_RNA.1.RNA.pdf"), width = 8, height = 16,
       plot = pheatmap(TotalCD8.LGAver$RNA[unique(RNAcelltypeheatgenes$gene),], scale = "row", border_color=NA, fontsize_row = 8))

#create Umap in the space of the selected genes
TotalCD8.LG <- RunUMAP(TotalCD8.LG, assay = "RNA", reduction = "pcaledge",
                    n.neighbors = 15, min.dist = 0.01, spread =  12, dims = 1:15)
TotalCD8.LGumap <- DimPlot(TotalCD8.LG, 
                        group.by = "RNA_snn_res.1", label = TRUE) + NoLegend()
TotalCD8.LGumap <- AugmentPlot(TotalCD8.LGumap, width = 4, height = 4, dpi = 300)
ggsave(TotalCD8.LGumap, filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD8.LGumap.WCTRNAclusters.pdf", width = 4, height = 4)

# calculate module scores of the pathways in the cells and plot
TotalCD8.LGpaths <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                      "HALLMARK_INFLAMMATORY_RESPONSE",
                      "btm_M4.1_cell cycle (I)",
                      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
                      "KEGG_RIBOSOME",
                      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                      "KEGG_OXIDATIVE_PHOSPHORYLATION",
                      "reactome_Fatty acid metabolism")
TotalCD8.LG <- AddModuleScore(TotalCD8.LG, features = genesets[TotalCD8.LGpaths], name = "pathway")
DotPlot(TotalCD8.LG, features =  c("pathway1","pathway2","pathway3","pathway4","pathway5")) + RotatedAxis()
ggsave(DotPlot(TotalCD8.LG, features =  c("pathway1","pathway2","pathway3","pathway4","pathway5","pathway6","pathway7","pathway8")) + RotatedAxis(),
       filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD8.DotPlot.pdf", width = 4, height = 4
)
FeaturePlot(TotalCD8.LG, features = "pathway3")
ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("JUN"),
                   order = TRUE, ncol = 3),
       filename = "TotalC8heatplot.JUN.WCTRNAclusters.pdf", width = 6, height = 6)
ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("MX1"),
                   order = TRUE, ncol = 3),
       filename = "TotalCD8.LGheatplot.MX1.WCTRNAclusters.pdf", width = 6, height = 6)
ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("MKI67"),
                   order = TRUE, ncol = 3),
       filename = "TotalCD8.LGheatplot.MKI67.WCTRNAclusters.pdf", width = 6, height = 6)
TotalCD8.LG <- AddModuleScore(TotalCD8.LG, features = cd8_exhaustion_list, name = "CD8_exhaustion")
ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("CD8_exhaustion1","CD8_exhaustion2","CD8_exhaustion3","CD8_exhaustion4","CD8_exhaustion5"),
                   order = TRUE, ncol = 3),
       filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalC8heatplot.Exhaustion.WCTRNAclusters.pdf", width = 18, height = 6)
ggsave(DotPlot(TotalCD8.LG, 
               features = c("CD8_exhaustion1","CD8_exhaustion2","CD8_exhaustion3","CD8_exhaustion4","CD8_exhaustion5")),
       filename = "withinWCTCelltypeRNAclustPathLeadingEdge/TotalC8Dotplot.Exhaustion.WCTRNAclusters.pdf", width = 18, height = 12)
TotalCD8.LG <- AddModuleScore(TotalCD8.LG, features = list(exhaustionscore = c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "TOX")), name = "exhaustionscore")
TotalCD8.LG <- ScaleData(TotalCD8.LG, features = c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "TOX"))
DoHeatmap(TotalCD8.LG,features = c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "TOX"))
pheatmap(TotalCD8.LGAver$RNA[c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "TOX"),])

ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("exhaustionscore1"),
                   order = TRUE),
       filename = "TotalC8heatplot.ExhaustionScoreWCTRNAclusters.pdf", width = 6, height = 6)
ggsave(FeaturePlot(TotalCD8.LG, 
                   features = c("PC1"),
                   order = TRUE),
       filename = "TotalC8heatplot.PC1.WCTRNAclusters.pdf", width = 10, height = 6)

ggsave(DotPlot(TotalCD8.LG, 
               features = c("CD8_exhaustion1","CD8_exhaustion2","CD8_exhaustion3","CD8_exhaustion4","CD8_exhaustion5","exhaustionscore1")),
       filename = "TotalC8Dotplot.Exhaustion.WCTRNAclusters.pdf", width = 18, height = 12)

TotalCD8.LGumapClass <- DimPlot(subset(TotalCD8.LG, Sorted == "N"), 
                             group.by = "WCTRNAclusters", label = TRUE, split.by = "Class") + NoLegend()
ggsave(TotalCD8.LGumapClass, filename = "TotalCD8.LGumapClass.WCTRNAclusters.pdf", width = 12, height = 6)
TotalCD8.LGumapSeverityOutcome <- DimPlot(subset(TotalCD8.LG, Sorted == "N"), 
                                       group.by = "WCTRNAclusters", label = TRUE, split.by = "severity.outcome") + NoLegend()
ggsave(TotalCD8.LGumapSeverityOutcome, filename = "TotalCD8.LGumapSeverity.WCTRNAcelltype.pdf", width = 12, height = 6)
FeaturePlot(subset(TotalCD8.LG, Sorted == "N"), features = c("MKI67"))
FeatureScatter(subset(TotalCD8.LG, Sorted =="N"), feature1 = "limmacite_CD25", "limmacite_CD127")
DotPlot(subset(TotalCD8.LG, Sorted == "N"), features = c("limmacite_CD25", "limmacite_CD127", "limmacite_HLA-DR",
                                                      "limmacite_CD71","limmacite_CD278","limmacite_CD69",
                                                      "limmacite_CD38", "limmacite_CD194", "limmacite_KLRG1",
                                                      "MKI67", "ISG15", "GZMA", "FOXP3", "TCF7", "TNFAIP3", "SELL", "CCR7"),
        group.by = "WCTRNAclusters") + RotatedAxis()

TotalCD8.LG.RNA.Labels = data.frame(WCTRNAclusters = TotalCD8.LG$WCTRNAclusters)
write.csv(TotalCD8.LG.RNA.Labels, file = "withinWCTCelltypeRNAclustPathLeadingEdge/RNAclustLabels/TotalCD8.LG.RNA.Labels.csv")
# TotalCD8.LG = AddModuleScore(TotalCD8.LG, features = genesets)

source("utilityFunctions/SubsetPlottingFunctions.r")
### make boxplots of the frequencies
TotalCD4.LG$WCTRNAclustParent <- rep("TotalCD4", length(colnames(TotalCD4.LG)))
TotalCD8.LG$WCTRNAclustParent <- rep("TotalCD8", length(colnames(TotalCD8.LG)))

remerge3.sng = merge(TotalCD4.LG, list(TotalCD8.LG))
remerge3.unsort = subset(remerge3.sng, subset = Sorted == "N")
remerge3.unsort.T0 = subset(remerge3.sng, subset = Timepoint == "T0")

# cell freq check of Seurat object
# get cell freq of both merged and course celltypes
remerge3.unsort_WCTRNAclusters <- data.frame(table(remerge3.unsort$sample_id, remerge3.unsort$WCTRNAclusters)) %>% 
  group_by(Var1) %>% mutate(ratio = Freq/sum(Freq))
remerge3.unsort_WCTparent <- data.frame(table(remerge3.unsort$sample_id, remerge3.unsort$WCTRNAclustParent)) %>% 
  group_by(Var1) %>% mutate(ratio = Freq/sum(Freq))

# add metadata
remerge3.unsort_WCTRNAclusters <- addcellmeta(remerge3.unsort_WCTRNAclusters)
remerge3.unsort_WCTparent <- addcellmeta(remerge3.unsort_WCTparent)
#set parents with less than 10 cells to NA
remerge3.unsort_WCTparent[which(remerge3.unsort_WCTparent$Freq < 10),c("Freq","ratio")] <- NA

# parse the dataframe to celltype matrix
library(data.table)
remerge3.unsort_WCTRNAclusters_mtx <- remerge3.unsort_WCTRNAclusters %>% 
  dcast(Var1+Timepoint+sample_id+Batch+severity+
          days_since_symptoms_onset+PC1class+PC1+Subject+severity_outcome ~ Var2, value.var = "ratio")

remerge3.unsort_WCTparent_mtx <- remerge3.unsort_WCTparent %>% 
  dcast(Var1 ~ Var2, value.var = "ratio")


# join merged and coruse together for getting freq of parent ratios
merge_cell_UnSort_WCTmerged_mtx <- left_join(remerge3.unsort_WCTRNAclusters_mtx, 
                                             remerge3.unsort_WCTparent_mtx, 
                                             "Var1")

merge_cell_UnSort_WCTmerged_to_parent_mtx <- merge_cell_UnSort_WCTmerged_mtx %>%
  mutate_at(vars(contains('TotalCD4.LG_')), ~(.)/TotalCD4) %>%
  mutate_at(vars(contains('TotalCD8.LG_')), ~(.)/TotalCD8)


merge_cell_UnSort_WCTmerged_to_parent <- merge_cell_UnSort_WCTmerged_to_parent_mtx %>% 
  melt(id = c("Var1", "Timepoint", "sample_id", "Batch", "severity", "days_since_symptoms_onset","PC1","PC1class","Subject","severity_outcome"))

library(plyr)
plot5cat(merge_cell_UnSort_WCTmerged_to_parent, width = 12, height = 8,
         "withinWCTCelltypeRNAclustPathLeadingEdge/merge.UnSort.WCTmerged.TotalCD4and8", "toparent")

saveRDS(merge_cell_UnSort_WCTmerged_to_parent, file = "TotalCD4and8.LG.20200824.merge_cell_UnSort_WCTmerged_to_parent.Rds")
saveRDS(merge_cell_UnSort_WCTmerged_to_parent_mtx, file = "TotalCD4and8.LG.20200824.merge_cell_UnSort_WCTmerged_to_parent_mtx.Rds")

write.csv(merge_cell_UnSort_WCTmerged_to_parent_mtx, file = "TotalCD4and8.LG.20200824.proportionToParent.csv")

### Create heatmaps of the CD4 and CD8 gene expression for the selected pathway genes
library(ComplexHeatmap)

library(viridis)
viridis(9)
dir.create("plots")

TotalCD4.LG = ScaleData(TotalCD4.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG,tO.PC1.SelectedPath.leadingEdge.LG)))
TotalCD4.selectclust.LG = subset(TotalCD4.LG, Sorted == "N" & RNA_snn_res.1 %in% c(14)& Timepoint == "T0") #

# create vector to label genesets
inflamgenes <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "HALLMARK_INFLAMMATORY_RESPONSE"
        
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "HALLMARK_INFLAMMATORY_RESPONSE"
          
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
    )

cellcyclegenes <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "btm_M4.1_cell cycle (I)"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
         
          "btm_M4.1_cell cycle (I)"
          
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

T1INFgenes <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "GO_RESPONSE_TO_TYPE_I_INTERFERON"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "GO_RESPONSE_TO_TYPE_I_INTERFERON"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

HALLMARK_TNFA_SIGNALING_VIA_NFKB <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

ribosome <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
       
        "KEGG_RIBOSOME"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_RIBOSOME"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

KEGG_GLYCOLYSIS_GLUCONEOGENESIS <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

KEGG_OXIDATIVE_PHOSPHORYLATION <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "KEGG_OXIDATIVE_PHOSPHORYLATION"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_OXIDATIVE_PHOSPHORYLATION"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

reactome_Fattyacidmetabolism <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "reactome_Fatty acid metabolism"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD4_Naive",
          "CD4_Mem",
          "Treg"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "reactome_Fatty acid metabolism"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD4_Naive",
            "CD4_Mem",
            "Treg"
          )  
      ),]$leadingEdge),  " ")))
  )

ha.CD4 = HeatmapAnnotation(
  # PC1 = TotalCD4.selectclust.LG$PC1,
  # PC1class = TotalCD4.selectclust.LG$PC1_cat,
  # Class = TotalCD4.selectclust.LG$Class,
  Cluster = TotalCD4.selectclust.LG$WCTRNAclusters
)

# list genes to be labeled in text on the plot
GOI = c("MX1","IFNAR1","IFNGR2","IRF3", "JUNB","PCNA","CEBPB","SOCS3","IL23A","EIF1","HIF1A","ALDOA", "POLA1",
        "STAT2","OAS1","NFKBIA","IRF7","MKI67", "TNFAIP3","TNF","TAP1","CD69","RELA", "SELL" ,"ISG15","AREG",
        "RPL13A","RPL22")
label <- which(rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% GOI)
label_genes <- rownames(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data"))[label]
  
  
library("circlize")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))

pdf("plots/CD4.Fig4Aheat.LG.geneslabelled.pdf", width = 10, height = 10)
Heatmap(as.matrix(GetAssayData(TotalCD4.selectclust.LG, assay = "RNA", slot = "scale.data")), name = "CD4", show_column_names = FALSE,
        top_annotation = ha.CD4, cluster_columns = FALSE, column_title_rot = 90,
        column_split = TotalCD4.selectclust.LG$WCTRNAclusters,
        col = col_fun,
        row_km = 8,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(inflamgenes + 0, name = "HALLMARK_INFLAMMATORY_RESPONSE", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(HALLMARK_TNFA_SIGNALING_VIA_NFKB + 0, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(cellcyclegenes + 0, name = "btm_M4.1_cellcycle", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(T1INFgenes + 0, name = "GO_RESPONSE_TO_TYPE_I_INTERFERON", col = c("0" = "white", "1" = "black"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(ribosome + 0, name = "KEGG_RIBOSOME", col = c("0" = "white", "1" = "orange"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_GLYCOLYSIS_GLUCONEOGENESIS + 0, name = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", col = c("0" = "white", "1" = "pink"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_OXIDATIVE_PHOSPHORYLATION + 0, name = "KEGG_OXIDATIVE_PHOSPHORYLATION", col = c("0" = "white", "1" = "green"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(reactome_Fattyacidmetabolism + 0, name = "reactome_Fatty acid metabolism", col = c("0" = "white", "1" = "grey"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 14)))
dev.off()

Idents(TotalCD4.LG) <- "RNA_snn_res.1"
TotalCD4.T0.LG = subset(TotalCD4.LG, Sorted == "N" & Timepoint == "T0")

TotalCD4.T0.aver.LG = AverageExpression(TotalCD4.T0.LG, return.seurat = TRUE)
TotalCD4.T0.aver.LG = ScaleData(TotalCD4.T0.aver.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG,tO.PC1.SelectedPath.leadingEdge.LG)))

GOI = c("MX1","IFNAR1","IFNGR2","IRF3", "JUNB","PCNA","CEBPB","SOCS3","IL23A","EIF1","HIF1A","ALDOA", "POLA1",
        "STAT2","OAS1","NFKBIA","IRF7","MKI67", "TNFAIP3","TNF","TAP1","CD69","RELA", "SELL" ,"ISG15","AREG","RPL13A","RPL22")
label <- which(rownames(GetAssayData(TotalCD4.T0.aver.LG, assay = "RNA", slot = "scale.data")) %in% GOI)
label_genes <- rownames(GetAssayData(TotalCD4.T0.aver.LG, assay = "RNA", slot = "scale.data"))[label]

pdf("plots/CD4.Fig4Aheat.AverLG.geneslabelled.pdf", width = 6.5, height = 6)
Heatmap(as.matrix(GetAssayData(TotalCD4.T0.aver.LG, assay = "RNA", slot = "scale.data")), name = "CD4", show_column_names = TRUE,
        cluster_columns = TRUE, column_title_rot = 90,
        #column_split = TotalCD8.T0.aver.LG$,
        col = col_fun,
        row_km = 8,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(inflamgenes + 0, name = "HALLMARK_INFLAMMATORY_RESPONSE", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(HALLMARK_TNFA_SIGNALING_VIA_NFKB + 0, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(cellcyclegenes + 0, name = "btm_M4.1_cellcycle", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(T1INFgenes + 0, name = "GO_RESPONSE_TO_TYPE_I_INTERFERON", col = c("0" = "white", "1" = "black"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(ribosome + 0, name = "KEGG_RIBOSOME", col = c("0" = "white", "1" = "orange"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_GLYCOLYSIS_GLUCONEOGENESIS + 0, name = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", col = c("0" = "white", "1" = "pink"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_OXIDATIVE_PHOSPHORYLATION + 0, name = "KEGG_OXIDATIVE_PHOSPHORYLATION", col = c("0" = "white", "1" = "green"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(reactome_Fattyacidmetabolism + 0, name = "reactome_Fatty acid metabolism", col = c("0" = "white", "1" = "grey"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 8)))
dev.off()

### heatmap for the CD8s
TotalCD8.LG = ScaleData(TotalCD8.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG.CD8,tO.PC1.SelectedPath.leadingEdge.LG.CD8)))
TotalCD8.selectclust.LG = subset(TotalCD8.LG, Sorted == "N" & RNA_snn_res.1 %in% c(14) & Timepoint == "T0")
# also label genesets
inflamgenes <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "HALLMARK_INFLAMMATORY_RESPONSE"
        
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "HALLMARK_INFLAMMATORY_RESPONSE"
          
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

cellcyclegenes <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "btm_M4.1_cell cycle (I)"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "btm_M4.1_cell cycle (I)"
          
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

T1INFgenes <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "GO_RESPONSE_TO_TYPE_I_INTERFERON"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "GO_RESPONSE_TO_TYPE_I_INTERFERON"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

HALLMARK_TNFA_SIGNALING_VIA_NFKB <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

ribosome <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "KEGG_RIBOSOME"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_RIBOSOME"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

KEGG_GLYCOLYSIS_GLUCONEOGENESIS <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

KEGG_OXIDATIVE_PHOSPHORYLATION <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "KEGG_OXIDATIVE_PHOSPHORYLATION"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "KEGG_OXIDATIVE_PHOSPHORYLATION"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

reactome_Fattyacidmetabolism <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% 
  union(unique(unlist(strsplit(as.character(t0_plus_t0_covid_only.PC1.fgseaatable[
    which(
      t0_plus_t0_covid_only.PC1.fgseaatable$pathway %in% c(
        
        "reactome_Fatty acid metabolism"
      ) 
      &
        t0_plus_t0_covid_only.PC1.fgseaatable$celltype %in% c(
          "CD8_Naive",
          "CD8_Mem"
        )  
    ),]$leadingEdge),  " ")))
    ,
    unique(unlist(strsplit(as.character(t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable[
      which(
        t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$pathway %in% c(
          
          "reactome_Fatty acid metabolism"
        ) 
        &
          t0_plus_healthy.healthy_vs_covid.COVID.Healthy.fgseaatable$celltype %in% c(
            "CD8_Naive",
            "CD8_Mem"
          )  
      ),]$leadingEdge),  " ")))
  )

ha.CD8 = HeatmapAnnotation(
  #PC1 = TotalCD8.selectclust.LG$PC1,
  #PC1class = TotalCD8.selectclust.LG$PC1_cat,
  #Class = TotalCD8.selectclust.LG$Class,
  Cluster = TotalCD8.selectclust.LG$WCTRNAclusters
)

GOI = c("MX1","IFNAR1","IFNGR2","IRF3", "JUNB","PCNA","CEBPB","SOCS3","IL23A","EIF1","HIF1A","ALDOA", "POLA1", "HELLS", "PRC1",
        "STAT2","OAS1","NFKBIA","IRF7","MKI67", "TNFAIP3","TNF","TAP1","CD69","RELA", "SELL" ,"ISG15","AREG","RPL13A","RPL22")
label <- which(rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")) %in% GOI)
label_genes <- rownames(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data"))[label]

library("circlize")
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#FF00FF", "#000000", "#FFFF00"))
pdf("plots/CD8.Fig4Aheat.geneslabelled.LG.pdf", width = 10, height = 10)
Heatmap(as.matrix(GetAssayData(TotalCD8.selectclust.LG, assay = "RNA", slot = "scale.data")), name = "CD8", show_column_names = FALSE,
        top_annotation = ha.CD8, cluster_columns = FALSE, column_title_rot = 90,
        column_split = TotalCD8.selectclust.LG$WCTRNAclusters,
        col = col_fun,
        row_km = 8,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(inflamgenes + 0, name = "HALLMARK_INFLAMMATORY_RESPONSE", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(HALLMARK_TNFA_SIGNALING_VIA_NFKB + 0, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(cellcyclegenes + 0, name = "btm_M4.1_cellcycle", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(T1INFgenes + 0, name = "GO_RESPONSE_TO_TYPE_I_INTERFERON", col = c("0" = "white", "1" = "grey"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(ribosome + 0, name = "KEGG_RIBOSOME", col = c("0" = "white", "1" = "orange"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_GLYCOLYSIS_GLUCONEOGENESIS + 0, name = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", col = c("0" = "white", "1" = "pink"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_OXIDATIVE_PHOSPHORYLATION + 0, name = "KEGG_OXIDATIVE_PHOSPHORYLATION", col = c("0" = "white", "1" = "green"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(reactome_Fattyacidmetabolism + 0, name = "reactome_Fatty acid metabolism", col = c("0" = "white", "1" = "black"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 14)))
dev.off()


Idents(TotalCD8.LG) <- "RNA_snn_res.1"
TotalCD8.T0.LG = subset(TotalCD8.LG, Sorted == "N" & Timepoint == "T0")
TotalCD8.T0.LG = subset(TotalCD8.T0.LG, idents = c("16","17"), invert = TRUE)
TotalCD8.T0.aver.LG = AverageExpression(TotalCD8.T0.LG, return.seurat = TRUE)
TotalCD8.T0.aver.LG = ScaleData(TotalCD8.T0.aver.LG, assay="RNA", features = unique(c(tO.COVIDvsHealthy.SelectedPath.leadingEdge.LG.CD8,tO.PC1.SelectedPath.leadingEdge.LG.CD8)))

GOI = c("MX1","IFNAR1","IFNGR2","IRF3", "JUNB","PCNA","CEBPB","SOCS3","IL23A","EIF1","HIF1A","ALDOA", "POLA1", "HELLS", "PRC1",
        "STAT2","OAS1","NFKBIA","IRF7","MKI67", "TNFAIP3","TNF","TAP1","CD69","RELA", "SELL" ,"ISG15","AREG","RPL13A","RPL22")
label <- which(rownames(GetAssayData(TotalCD8.T0.aver.LG, assay = "RNA", slot = "scale.data")) %in% GOI)
label_genes <- rownames(GetAssayData(TotalCD8.T0.aver.LG, assay = "RNA", slot = "scale.data"))[label]

pdf("plots/CD8.Fig4Aheat.geneslabelled.AverLG.pdf", width = 6.5, height = 6)
Heatmap(as.matrix(GetAssayData(TotalCD8.T0.aver.LG, assay = "RNA", slot = "scale.data")), name = "CD8", show_column_names = TRUE,
        cluster_columns = TRUE, column_title_rot = 90,
        #column_split = TotalCD8.T0.aver.LG$,
        col = col_fun,
        row_km = 8,
        row_names_gp = gpar(fontsize = 10)
) +
  Heatmap(inflamgenes + 0, name = "HALLMARK_INFLAMMATORY_RESPONSE", col = c("0" = "white", "1" = "purple"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(HALLMARK_TNFA_SIGNALING_VIA_NFKB + 0, name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", col = c("0" = "white", "1" = "red"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(cellcyclegenes + 0, name = "btm_M4.1_cellcycle", col = c("0" = "white", "1" = "blue"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(T1INFgenes + 0, name = "GO_RESPONSE_TO_TYPE_I_INTERFERON", col = c("0" = "white", "1" = "black"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(ribosome + 0, name = "KEGG_RIBOSOME", col = c("0" = "white", "1" = "orange"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_GLYCOLYSIS_GLUCONEOGENESIS + 0, name = "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", col = c("0" = "white", "1" = "pink"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(KEGG_OXIDATIVE_PHOSPHORYLATION + 0, name = "KEGG_OXIDATIVE_PHOSPHORYLATION", col = c("0" = "white", "1" = "green"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  Heatmap(reactome_Fattyacidmetabolism + 0, name = "reactome_Fatty acid metabolism", col = c("0" = "white", "1" = "grey"),
          show_heatmap_legend = FALSE, width = unit(2, "mm"), column_names_gp = gpar(fontsize = 6)) +
  rowAnnotation(foo = anno_mark(at = label,
                              labels = label_genes,
                              labels_gp = gpar(fontsize = 8)))
dev.off()

pdf("plots/CD8.Fig4.Vlnplot.pdf", width = 10, height = 5)
VlnPlot(TotalCD8.T0.LG, features = c("limmacite_CD38", "limmacite_HLA-DR","limmacite_CD278","limmacite_CD86"), group.by = "RNA_snn_res.1", pt.size = 0, ncol = 4)
dev.off()

TotalCD4.T0.LG$RNA_snn_res.1.num = factor(Idents(TotalCD4.T0.LG), levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
pdf("plots/CD4.Fig4.Vlnplot.pdf", width = 10, height = 5)
VlnPlot(subset(TotalCD4.T0.LG, RNA_snn_res.1 == "<NA>", invert=TRUE), features = c("limmacite_CD38", "limmacite_HLA-DR","limmacite_CD278","limmacite_CD71"), group.by = "RNA_snn_res.1.num", pt.size = 0, ncol = 4)
dev.off()

TotalCD4.T0.LGAver = AverageExpression_MeanOnly(TotalCD4.T0.LG, return.seurat=F)
mat_breaks.temp <- quantile_breaks(as.matrix(TotalCD4.T0.LGAver$limmaCITE[c("CD38", "HLA-DR","CD278","CD279", "CD71","KLRG1", "CD28","CD45RA","CD45RO"),]), n = 101)
ggsave(filename = paste0("plots/","TotalCD4.T0.LG.Aver.SelectProts.pdf"), width = 8, height = 4,
       plot = pheatmap(TotalCD4.T0.LGAver$limmaCITE[c("CD38", "HLA-DR","CD278","CD279", "CD71","KLRG1", "CD28","CD45RA","CD45RO"),], scale = "none", border_color=NA, viridis(n=100), breaks = mat_breaks.temp))

TotalCD8.T0.LGAver = AverageExpression_MeanOnly(TotalCD8.T0.LG, return.seurat=F)
mat_breaks.temp <- quantile_breaks(as.matrix(TotalCD8.T0.LGAver$limmaCITE[c("CD38", "HLA-DR","CD278","CD279", "CD86","KLRG1", "CD28","CD45RA","CD45RO"),]), n = 101)
ggsave(filename = paste0("plots/","TotalCD8.T0.LG.Aver.SelectProts.pdf"), width = 8, height = 4,
       plot = pheatmap(TotalCD8.T0.LGAver$limmaCITE[c("CD38", "HLA-DR","CD278","CD279", "CD86","KLRG1", "CD28","CD45RA","CD45RO"),], scale = "none", border_color=NA, viridis(n=100), breaks = mat_breaks.temp))

saveRDS(TotalCD8.LG, "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD8.LG.Rds")
saveRDS(TotalCD4.LG, "withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LG.Rds")
TotalCD8.LG = readRDS("withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD8.LG.Rds")
TotalCD4.LG = readRDS("withinWCTCelltypeRNAclustPathLeadingEdge/TotalCD4.LG.Rds")

CD4.T0.LG.SpecificGatingPerCluster = prop.table(table(Idents(TotalCD4.T0.LG), TotalCD4.T0.LG$specific_gating), margin = 1)
CD4.T0.LG.GateGroupMergedPerCluster = prop.table(table(Idents(TotalCD4.T0.LG), TotalCD4.T0.LG$gate_group_merged), margin = 1)
CD4.T0.LG.WCTmergedPerCluster = prop.table(table(Idents(TotalCD4.T0.LG), TotalCD4.T0.LG$WCTmergedcelltype), margin = 1)
CD4.T0.LG.WCTmergedPerCluster = prop.table(table(Idents(TotalCD4.T0.LG), TotalCD4.T0.LG$WCTmergedcelltype), margin = 1)

colnames(CD4.T0.LG.WCTmergedPerCluster)[17] <- "Treg.ActT"
colnames(CD4.T0.LG.WCTmergedPerCluster)
colnames(CD4.T0.LG.SpecificGatingPerCluster)
colnames(CD4.T0.LG.GateGroupMergedPerCluster)
pheatmap(CD4.T0.LG.WCTmergedPerCluster[,which(!colnames(CD4.T0.LG.WCTmergedPerCluster) %in% 
                                                   c("CD4_Mem_CD69pos","CD4_Mem_Activated.HLADRhi",
                                                     "CD4_Mem_CD22hi","CD4_Mem_MAIT",
                                                     "CD4_Naive_MAIT","CD4_Naive_CD38hi",
                                                     "CD4_Mem_integrin.hi","CD4_Mem_CD41hi",
                                                     "CD4_Mem_IgG.IsoBinding","CD4_Mem_CD146hi",
                                                     "CD4_Naive_IgG.IsoBinding"))],
         border_color = NA, filename = "plots/OlapCD4T0.LGclusts.pdf", width = 3, height = 6)
pheatmap(cbind(Tfh = CD4.T0.LG.SpecificGatingPerCluster[,2], GatedTreg = CD4.T0.LG.GateGroupMergedPerCluster[,13]),
         border_color = NA, filename = "plots/OlapCD4T0.LGclusts.GatedPop.pdf", width = 2, height = 6)


CD8.T0.LG.SpecificGatingPerCluster = prop.table(table(Idents(TotalCD8.T0.LG), TotalCD8.T0.LG$specific_gating), margin = 1)
CD8.T0.LG.GateGroupMergedPerCluster = prop.table(table(Idents(TotalCD8.T0.LG), TotalCD8.T0.LG$gate_group_merged), margin = 1)
CD8.T0.LG.WCTmergedPerCluster = prop.table(table(Idents(TotalCD8.T0.LG), TotalCD8.T0.LG$WCTmergedcelltype), margin = 1)
colnames(CD8.T0.LG.WCTmergedPerCluster)
pheatmap(CD8.T0.LG.WCTmergedPerCluster[,which(!colnames(CD8.T0.LG.WCTmergedPerCluster) %in% 
                                             c("CD8_Mem_CD158hi","CD8_Naive_CD41hi",
                                               "CD8_Mem_CD41hi","CD8_Mem_CD69pos",
                                               "CD8_Naive_MAIT","CD8_Mem_IgG.IsoBinding",
                                               "CD4_Naive","CD8_Naive_NKT","CD8_Naive_CD101hi",
                                               "CD8_Mem_NKT","CD8_Mem_CD101hi"))],
         border_color = NA, filename = "plots/OlapCD8T0.LGclusts.pdf", width = 4, height = 6)


#import list of expanded T cells
expandedT = read.table("Metadata/2020_07_28_tcell_expansion_df.tsv", header = TRUE)

TotalCD4.LG$Expanded <- ifelse(colnames(TotalCD4.LG) %in% 
                                 expandedT$barcodeBatch[which(expandedT$expanded == TRUE)], TRUE, FALSE)
prop.table(table(TotalCD4.LG$RNA_snn_res.1, TotalCD4.LG$Expanded), margin = 1)

TotalCD4.T0.LG$Expanded <- ifelse(colnames(TotalCD4.T0.LG) %in% 
                                 expandedT$barcodeBatch[which(expandedT$expanded == TRUE)], TRUE, FALSE)
prop.table(table(TotalCD4.T0.LG$RNA_snn_res.1.num, TotalCD4.T0.LG$Expanded), margin = 1)

TotalCD8.LG$Expanded <- ifelse(colnames(TotalCD8.LG) %in% 
                                    expandedT$barcodeBatch[which(expandedT$expanded == TRUE)], TRUE, FALSE)

TotalCD8.T0.LG$Expanded <- ifelse(colnames(TotalCD8.T0.LG) %in% 
                                    expandedT$barcodeBatch[which(expandedT$expanded == TRUE)], TRUE, FALSE)
prop.table(table(TotalCD8.T0.LG$RNA_snn_res.1, TotalCD8.T0.LG$Expanded), margin = 1)

TotalCD4.T0.LG.unsort <- subset(TotalCD4.T0.LG, Sorted == "N")
TotalCD8.T0.LG.unsort <- subset(TotalCD8.T0.LG, Sorted == "N")

#make barplots of expanded cell percentage per cluster
pdf(file = "plots/CD4.ExpandedPercluster.pdf", width = 4, height = 4)
plot(as.numeric(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1.num)), 
     as.numeric(prop.table(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1.num, TotalCD4.T0.LG.unsort$Expanded), margin = 1)[,2]),
     xlab = "Total Cells in Cluster", ylab = "Proportion of Cluster Expanded", main = "CD4")
text(as.numeric(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1.num)), 
     as.numeric(prop.table(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1.num, TotalCD4.T0.LG.unsort$Expanded), margin = 1)[,2]),
     labels = names(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1.num)), pos = 4)
dev.off()

pdf(file = "plots/CD8.ExpandedPercluster.pdf", width = 4, height = 4)
plot(as.numeric(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1)), 
     as.numeric(prop.table(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1, TotalCD8.T0.LG.unsort$Expanded), margin = 1)[,2]),
     xlab = "Total Cells in Cluster", ylab = "Proportion of Cluster Expanded", main = "CD8")
text(as.numeric(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1)), 
     as.numeric(prop.table(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1, TotalCD8.T0.LG.unsort$Expanded), margin = 1)[,2]),
     labels = names(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1)), pos = 4)
dev.off()

table(subset(TotalCD4.T0.LG.unsort, RNA_snn_res.1 == "7")$Class, subset(TotalCD4.T0.LG.unsort, RNA_snn_res.1 == "7")$Expanded)

pdf(file = "plots/CD4.ExpandedPercluster.Barplot.pdf", width = 8, height = 4)
barplot(prop.table(table(TotalCD4.T0.LG.unsort$RNA_snn_res.1, TotalCD4.T0.LG.unsort$Expanded), 
                   margin = 1)[c(7,10,12,3,13,8,2,1,11,4,10,5,14,15,6),2],
        horiz = FALSE,
        ylim = c(0,1)
        )
dev.off()
pdf(file = "plots/CD8.ExpandedPercluster.Barplot.pdf", width = 8, height = 4)
barplot(prop.table(table(TotalCD8.T0.LG.unsort$RNA_snn_res.1, TotalCD8.T0.LG.unsort$Expanded), 
                   margin = 1)[c(15,12,14,11,1,5,4,2,8,13,7,3,10,16,6,9),2],
        horiz = FALSE,
        ylim = c(0,1)
        )
dev.off()

prop.table(table(subset(TotalCD4.LG, idents = 7)$Class, subset(TotalCD4.LG, idents = 7)$Expanded),margin =1)
prop.table(table(subset(TotalCD4.LG, idents = 14)$Class, subset(TotalCD4.LG, idents = 14)$Expanded),margin =1)
prop.table(table(subset(TotalCD8.LG, idents = 6)$Class, subset(TotalCD8.LG, idents = 6)$Expanded),margin =1)
prop.table(table(subset(TotalCD8.LG, idents = 14)$Class, subset(TotalCD8.LG, idents = 14)$Expanded),margin =1)

# save objects
saveRDS(TotalCD4.LG, "TotalCD4.LG.20200914.Rds")
saveRDS(TotalCD4.T0.LG, "TotalCD4.T0.LG.20200914.Rds")
saveRDS(TotalCD8.LG, "TotalCD8.LG.20200914.Rds")
saveRDS(TotalCD8.T0.LG, "TotalCD8.T0.LG.20200914.Rds")




