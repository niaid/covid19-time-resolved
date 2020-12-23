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

# load data from Schulte-Schrepping et al. Cell. 2020 Sep 17; 182(6): 1419–1440.e23 PMID: 32810438
# Data was obtained from: https://www.fastgenomics.org
SchultePBMC <- readRDS("externalData/SchulteSchrepping/seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds")
colnames(SchultePBMC@meta.data)

merge.SNG <- readRDS("SeuratObjects/AllBatches_SeuratObj.rds")

pbmc.list = list(merge.SNG = merge.SNG, SchultePBMC = SchultePBMC)

for (i in names(pbmc.list)) {
  pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
}
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)

reference_dataset <- which(names(pbmc.list) == "merge.SNG")

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "LogNormalize", 
                                       anchor.features = pbmc.features, reference = reference_dataset,
                                       reduction = 'cca', dims = 1:30, assay = c("RNA","RNA"))
saveRDS(pbmc.anchors, "SeuratObjects/Schulte.pbmc.anchors.Rds")
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "LogNormalize")
DefaultAssay(pbmc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)

pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)

p1 <- DimPlot(pbmc.integrated, reduction = "umap", group.by = "id.celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p2 <- DimPlot(pbmc.integrated, reduction = "umap", group.by = "WCTcoursecelltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2


pbmcT.anchors <- FindTransferAnchors(reference = merge.SNG, query = SchultePBMC, 
                                        dims = 1:30, normalization.method = "LogNormalize",
                                      reduction = "cca" )
predictions <- TransferData(anchorset = pbmcT.anchors, refdata = merge.SNG$WCTcoursecelltype, 
                            dims = 1:30, weight.reduction = "cca")
SchultePBMC <- AddMetaData(SchultePBMC, metadata = predictions)
pheatmap(table(SchultePBMC$id.celltype, SchultePBMC$predicted.id), display_numbers = FALSE,
         filename = "SchultePBMC.predicted.idOlaps.pdf", width = 10, height = 10, scale = "row")
pheatmap(t(table(SchultePBMC$id.celltype, SchultePBMC$predicted.id)), display_numbers = FALSE,
           filename = "SchultePBMC.predicted.idOlaps.Transpose.pdf", width = 10, height = 10, scale = "row")
saveRDS(predictions, "SeuratObjects/predictedLabels.Schulte.Rds")
saveRDS(pbmcT.anchors, "SeuratObjects/Schulte.pbmc.anchors.Rds")


SchulteCohort2 <- readRDS("externalData/SchulteSchrepping/seurat_COVID19_PBMC_jonas_FG_2020-07-23.rds")

pbmcT.anchors.C2 <- FindTransferAnchors(reference = merge.SNG, query = SchulteCohort2, 
                                     dims = 1:30, normalization.method = "LogNormalize",
                                     reduction = "cca" )
predictions.C2 <- TransferData(anchorset = pbmcT.anchors.C2, refdata = merge.SNG$WCTcoursecelltype, 
                            dims = 1:30, weight.reduction = "cca")
SchulteCohort2 <- AddMetaData(SchulteCohort2, metadata = predictions.C2)
pheatmap(table(SchulteCohort2$cluster_labels_res.0.4, SchulteCohort2$predicted.id), display_numbers = FALSE,
         filename = "SchulteCohort2.predicted.idOlaps.pdf", width = 10, height = 10, scale = "row")
pheatmap(t(table(SchulteCohort2$cluster_labels_res.0.4, SchulteCohort2$predicted.id)), display_numbers = FALSE,
         filename = "SchulteCohort2.predicted.idOlaps.Transpose.pdf", width = 10, height = 10, scale = "row")

saveRDS(predictions.C2, "SeuratObjects/predictedLabels.Schulte.C2.Rds")
saveRDS(pbmcT.anchors.C2, "SeuratObjects/pbmcT.anchors.C2.Rds")
