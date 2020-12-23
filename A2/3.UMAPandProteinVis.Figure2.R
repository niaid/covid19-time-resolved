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

####UMAP plots for Figure 2
## import the seurat object downloaded from GEO
merge.SNG = readRDS("SeuratObjects/AllBatches_SeuratObj.rds")

##remove the CHI014 technical control sample
merge.SNG <- subset(merge.SNG, Donor == "CHI014", invert = TRUE)

## make a object of the cluster averages for heatmaps 
Idents(merge.SNG) <- "WCTcoursecelltype"
source("utilityFunctions/AverageExpression_MeanOnly.r")
environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions
merge.SNG.Aver = AverageExpression_MeanOnly(merge.SNG, assays = "limmaCITE")
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

f2 <- pOverA(0.01, 3)
ffun2 <- filterfun(f2)
merge.SNGADTfilt = genefilter(GetAssayData(merge.SNG[["limmaCITE"]], slot = "data"), ffun2)

mat_breaks <- quantile_breaks(as.matrix(merge.SNG.Aver$limmaCITE[which(merge.SNGADTfilt),]), n = 101)
pheatmap(merge.SNG.Aver$limmaCITE[which(merge.SNGADTfilt),which(!(colnames(merge.SNG.Aver$limmaCITE) %in% c("Dblt","Unknown","gated_out","dim","Tcell","dblt")))], scale = "none", border_color=NA, viridis(n=100), breaks = mat_breaks,
                       file = "Plots/COVID.WCTclust_CITE.aver.pdf", width = 6, height = 15)

innate.merge.SNG <- subset(merge.SNG, subset = WCTcoursecelltype %in% 
                             c("Mono_Classical", "Mono_NonClassical", "NK_CD16hi",
                               "NK_CD56loCD16lo", "pDC", "cDC", "Platelets", "NK_CD56hiCD16lo",
                               "Mono_Intermediate", "Granulocytes"))
adaptive.merge.SNG <- subset(merge.SNG, subset = WCTcoursecelltype %in% 
                             c("B_Naive", "B_Mem", "gammadeltaT",
                               "PB_Plasmablasts", "Treg", "CD4_Mem", "CD8_Mem", "CD8_Naive",
                               "CD4_Naive", "MAIT", "TissueResMemT", "DPT", "DNT", "TCRVbeta13.1pos"))

innate.merge.SNGADTfilt = genefilter(GetAssayData(innate.merge.SNG[["limmaCITE"]], slot = "data"), ffun2)
innate.merge.SNG <- RunUMAP(innate.merge.SNG, assay = "limmaCITE", features = rownames(innate.merge.SNG[["limmaCITE"]][which(innate.merge.SNGADTfilt),]),
                     n.neighbors = 15, min.dist = 0.01, spread =  5)
innate.p <- DimPlot(innate.merge.SNG, 
             group.by = "WCTmergedcelltype", label = FALSE)
innate.p.aug <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTmergedcelltype.labels.pdf", width = 8, height = 6)
ggsave(innate.p.aug, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTmergedcelltype.pdf", width = 6, height = 6)
innate.p <- DimPlot(innate.merge.SNG, 
                      group.by = "WCTcoursecelltype", label = TRUE) + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.pdf", width = 6, height = 6)
innate.p <- DimPlot(innate.merge.SNG, 
                    group.by = "WCTcoursecelltype", label = FALSE) + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.unlabelled.pdf", width = 6, height = 6)
innate.p <- FeaturePlot(innate.merge.SNG, 
                    features = "days_of_symptoms_onset") + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.days_of_symptoms_onset.pdf", width = 6, height = 6)
innate.p <- FeaturePlot(innate.merge.SNG, 
                        features = "cite_CD163") + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.CD163.pdf", width = 6, height = 6)
innate.p <- FeaturePlot(innate.merge.SNG, 
                        features = "cite_CD163") + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.CD163.pdf", width = 6, height = 6)
innate.p <- FeaturePlot(innate.merge.SNG, 
                        features = "cite_HLA-DR") + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "UMAP.innate.merge.SNG.filt.WCTcoursecelltype.HLADR.pdf", width = 6, height = 6)
innate.p <- FeaturePlot(innate.merge.SNG, 
                        features = "cite_CD184") + NoLegend()
innate.p <- AugmentPlot(innate.p, width = 6, height = 6)
ggsave(innate.p, filename = "Plots/UMAP.innate.merge.SNG.filt.WCTcoursecelltype.CD184.pdf", width = 6, height = 6)

adaptive.merge.SNGADTfilt = genefilter(GetAssayData(adaptive.merge.SNG[["limmaCITE"]], slot = "data"), ffun2)
adaptive.merge.SNG <- RunUMAP(adaptive.merge.SNG, assay = "limmaCITE", features = rownames(adaptive.merge.SNG[["limmaCITE"]][which(adaptive.merge.SNGADTfilt),]),
                            n.neighbors = 15, min.dist = 0.01, spread =  5)
adaptive.p <- DimPlot(adaptive.merge.SNG, 
                    group.by = "WCTmergedcelltype", label = FALSE)
adaptive.p.aug <- AugmentPlot(adaptive.p, width = 6, height = 6)
ggsave(adaptive.p, filename = "UMAP.adaptive.merge.SNG.filt.WCTmergedcelltype.labels.pdf", width = 8, height = 6)
ggsave(adaptive.p.aug, filename = "Plots/UMAP.adaptive.merge.SNG.filt.WCTmergedcelltype.pdf", width = 6, height = 6)
adaptive.p <- DimPlot(adaptive.merge.SNG, 
                      group.by = "WCTcoursecelltype", label = TRUE) + NoLegend()
adaptive.p <- AugmentPlot(adaptive.p, width = 6, height = 6)
ggsave(adaptive.p, filename = "Plots/UMAP.adaptive.merge.SNG.filt.WCTcoursecelltype.pdf", width = 6, height = 6)
adaptive.p <- DimPlot(adaptive.merge.SNG, 
                      group.by = "WCTcoursecelltype", label = FALSE) + NoLegend()
adaptive.p <- AugmentPlot(adaptive.p, width = 6, height = 6)
ggsave(adaptive.p, filename = "Plots/UMAP.adaptive.merge.SNG.filt.WCTcoursecelltype.unlabelled.pdf", width = 6, height = 6)
adaptive.p <- FeaturePlot(adaptive.merge.SNG, 
                      features = "days_of_symptoms_onset") + NoLegend()
adaptive.p <- AugmentPlot(adaptive.p, width = 6, height = 6)
ggsave(adaptive.p, filename = "Plots/UMAP.adaptive.merge.SNG.filt.WCTcoursecelltype.days_of_symptoms_onset.pdf", width = 6, height = 6)

Idents(innate.merge.SNG) <- "WCTmergedcelltype"
innate.merge.SNG.Aver = AverageExpression_MeanOnly(innate.merge.SNG, assays = "limmaCITE")
innate.mat_breaks <- quantile_breaks(as.matrix(innate.merge.SNG.Aver$limmaCITE[which(innate.merge.SNGADTfilt),c(1:11,13:20,22:25)]), n = 101)
pheatmap(innate.merge.SNG.Aver$limmaCITE[which(innate.merge.SNGADTfilt),c(1:11,13:20,22:25)], scale = "none", border_color=NA, viridis(n=100), breaks = innate.mat_breaks,
         file = "Plots/COVID.innate.WCTclust_CITE.aver.pdf", width = 6, height = 15)

Idents(adaptive.merge.SNG) <- "WCTmergedcelltype"
adaptive.merge.SNG.Aver = AverageExpression_MeanOnly(adaptive.merge.SNG, assays = "limmaCITE")
adaptive.mat_breaks <- quantile_breaks(as.matrix(adaptive.merge.SNG.Aver$limmaCITE[which(adaptive.merge.SNGADTfilt),c(1:69,71:73)]), n = 101)
pheatmap(adaptive.merge.SNG.Aver$limmaCITE[which(adaptive.merge.SNGADTfilt),c(1:69,71:73)], scale = "none", border_color=NA, viridis(n=100), breaks = adaptive.mat_breaks,
         file = "Plots/COVID.adaptive.WCTclust_CITE.aver.pdf", width = 12, height = 15)

Mono_Classical <- subset(merge.SNG, subset = WCTcoursecelltype %in% c( "Mono_Classical"))
Mono_Classical <- FindVariableFeatures(Mono_Classical, assay = "RNA")
Mono_Classical <- ScaleData(Mono_Classical)
Mono_Classical <- RunPCA(Mono_Classical)
ElbowPlot(Mono_Classical)
Mono_Classical <- RunUMAP(Mono_Classical, dims = 1:15)

MonoUMAP.Donor <- DimPlot(Mono_Classical, group.by = "Donor") + NoLegend()
ggsave(MonoUMAP.Donor, filename = "Plots/MonoUMAP.Donor.pdf", width = 6, height = 6)
MonoUMAP.Donor <- AugmentPlot(MonoUMAP.Donor, width = 6, height = 6)
ggsave(MonoUMAP.Donor, filename = "Plots/MonoUMAP.Donor.Aug.pdf", width = 6, height = 6)

MonoUMAP.Class <- DimPlot(Mono_Classical, group.by = "Class", cols = c("#00BA38", "#F8766D"))
ggsave(MonoUMAP.Class, filename = "Plots/MonoUMAP.Class.pdf", width = 7, height = 6)
MonoUMAP.Class <- AugmentPlot(MonoUMAP.Class, width = 6, height = 6)
ggsave(MonoUMAP.Class, filename = "Plots/MonoUMAP.Class.Aug.pdf", width = 6, height = 6)

MonoUMAP.DSM_cat <- DimPlot(subset(Mono_Classical, PC1_cat %in% c("PC1_low","PC1_high")), group.by = "PC1_cat",cols = c("#619CFF", "#F8766D"))
ggsave(MonoUMAP.DSM_cat, filename = "Plots/MonoUMAP.DSM_cat.pdf", width = 7, height = 6)
MonoUMAP.DSM_cat <- AugmentPlot(MonoUMAP.DSM_cat, width = 6, height = 6)
ggsave(MonoUMAP.DSM_cat, filename = "Plots/MonoUMAP.DSM_cat.Aug.pdf", width = 6, height = 6)

MonoUMAP.dayssinceonset <- FeaturePlot(subset(Mono_Classical, PC1_cat %in% c("PC1_low","PC1_high") & days_since_onset != "NA"), features = "days_since_onset")
ggsave(MonoUMAP.dayssinceonset, filename = "Plots/MonoUMAP.dayssinceonset.pdf", width = 7, height = 6)
MonoUMAP.dayssinceonset <- AugmentPlot(MonoUMAP.dayssinceonset, width = 6, height = 6)
ggsave(MonoUMAP.dayssinceonset, filename = "Plots/MonoUMAP.dayssinceonset.Aug.pdf", width = 6, height = 6)

MonoUMAP.CD163 <- DimPlot(subset(Mono_Classical, WCTmergedcelltype %in% c("Mono_Classical","Mono_Classical_CD163hi")), group.by = "WCTmergedcelltype", cols = c("grey","black"))
ggsave(MonoUMAP.CD163, filename = "Plots/MonoUMAP.CD163.pdf", width = 8, height = 6)
MonoUMAP.CD163 <- AugmentPlot(MonoUMAP.CD163, width = 6, height = 6)
ggsave(MonoUMAP.CD163, filename = "Plots/MonoUMAP.CD163.Aug.pdf", width = 6, height = 6)

MonoUMAP.CD163x <- FeaturePlot(subset(Mono_Classical, WCTmergedcelltype %in% c("Mono_Classical","Mono_Classical_CD163hi")), features = "cite_CD163",max.cutoff = 7.5 )
ggsave(MonoUMAP.CD163x, filename = "Plots/MonoUMAP.CD163expression.pdf", width = 6, height = 6)
MonoUMAP.CD163x <- AugmentPlot(MonoUMAP.CD163x, width = 6, height = 6)
ggsave(MonoUMAP.CD163x, filename = "Plots/MonoUMAP.CD163expression.Aug.pdf", width = 6, height = 6)

MonoUMAP.HLADRx <- FeaturePlot(subset(Mono_Classical, WCTmergedcelltype %in% c("Mono_Classical","Mono_Classical_CD163hi")), features = "cite_HLA-DR",max.cutoff = 7.5)
ggsave(MonoUMAP.HLADRx, filename = "Plots/MonoUMAP.HLA-DRexpression.pdf", width = 6, height = 6)
MonoUMAP.HLADRx <- AugmentPlot(MonoUMAP.HLADRx, width = 6, height = 6)
ggsave(MonoUMAP.HLADRx, filename = "Plots/MonoUMAP.HLA-DRexpression.Aug.pdf", width = 6, height = 6)
