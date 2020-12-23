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

#breaks function for use in plots
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
source("utilityFunctions/AverageExpression_MeanOnly.r")
environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions

######  Clustering based on protein (CITE) data
merge.SNG = readRDS(file = "SeuratObjects/merge.SNG.forClust.rds")

batches = unique(merge.SNG$Batch)
# make a list of all the objects per batch 
batchObjects = list()
for(i in 1:length(batches)){
  batchObjects[[i]] = subset(merge.SNG, subset = Batch %in% batches[i])
}

names(batchObjects) <- batches

# cluster within each batch
for(i in 1:length(batches)){
  CITEtemp = GetAssayData(batchObjects[[i]][["CITE"]], slot = "data")[c(1:79,84:192),]
  adt.dist.temp <- parDist(t(CITEtemp), threads = 16)
  batchObjects[[i]][["adt_snn"]] <- FindNeighbors(adt.dist.temp, nn.eps = 1)$snn
  batchObjects[[i]] <- FindClusters(batchObjects[[i]], resolution = c(0.5), graph.name = "adt_snn", algorithm = 1)
  batchObjects[[i]] <- DietSeurat(batchObjects[[i]])
}
##relabel BatchClusters
for(i in 1:length(batches)){
  batchObjects[[i]]$BatchClusters <- paste0(batchObjects[[i]]$Batch, "_", batchObjects[[i]]$adt_snn_res.0.5)
}
merge.SNG <- merge(batchObjects[[1]], batchObjects[2:length(batches)])
Idents(merge.SNG) <- "BatchClusters"
BatchClusterAverages = AverageExpression_MeanOnly(merge.SNG, return.seurat = TRUE)

## plot averages for cluter naming
mat_breaks.temp <- quantile_breaks(as.matrix(GetAssayData(BatchClusterAverages[["CITE"]])), n = 101)
ggsave(filename = "Plots/BatchClusterAverages.COVID_CITE.pdf", width = 18, height = 24,
         plot = pheatmap(GetAssayData(BatchClusterAverages[["CITE"]]), scale = "none", border_color=NA, 
                         inferno(length(mat_breaks.temp) - 1), breaks = mat_breaks.temp))


#name the clusters. bring back in the clusternames here from csv files
withinBatchClusterMergedcelltypeNames = read.csv("Metadata/withinBatchClusterMergedcelltypeNames.csv", header=TRUE)
Idents(merge.SNG) <- "BatchClusters"
cluster.ids <- withinBatchClusterMergedcelltypeNames$BatchClusters
setdiff(Idents(merge.SNG), cluster.ids)
setdiff(cluster.ids, Idents(merge.SNG))
merge.SNG$mergedcelltype <- plyr::mapvalues(x = Idents(merge.SNG), from = cluster.ids, to = as.character(withinBatchClusterMergedcelltypeNames$mergedcelltype))
merge.SNG$coursecelltype <- plyr::mapvalues(x = Idents(merge.SNG), from = cluster.ids, to = as.character(withinBatchClusterMergedcelltypeNames$coursecelltype))

#add metadata
samplemetadata = read.csv("Metadata/allbatches.HTOandTimepointMetadata.csv", header = TRUE, colClasses = c(rep("character", 7)))
merge.SNG$Batch_Sample = merge.SNG$Sample
merge.SNG$Sample = paste0(merge.SNG$Subject, "_", merge.SNG$autoHashcalls)
samplenames = unique(merge.SNG$Sample)
setdiff(samplemetadata$Sample, samplenames)
setdiff(samplenames, samplemetadata$Sample)
Idents(merge.SNG) <- "Sample"
merge.SNG$Age <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Age)
merge.SNG$Gender <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Gender)
merge.SNG$Ward <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Ward)
merge.SNG$Status <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Status)
merge.SNG$Pool <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Pool)
merge.SNG$Timepoint <- plyr::mapvalues(x = Idents(merge.SNG), from = samplemetadata$x, to = samplemetadata$Timepoint)
table(merge.SNG$x, merge.SNG$Pool) 
table(merge.SNG$x, merge.SNG$Age) 
table(merge.SNG$x, merge.SNG$Gender) 
table(merge.SNG$autoHashcalls, merge.SNG$Donor)
merge.SNG <- subset(merge.SNG, subset = Pool == "misassigned", invert = TRUE)

####re-assign some cells based on gating
table(merge.SNG$Batch, merge.SNG$mergedcelltype)
Idents(merge.SNG) <- "coursecelltype"
adt <- as.matrix(merge.SNG@assays$CITE@data)
merge.SNG$adjustedcelltype <- merge.SNG$mergedcelltype
merge.SNG$Batch <- as.character(merge.SNG$Batch)
merge.SNG$coursecelltype <- as.character(merge.SNG$coursecelltype)
adjustedcelltype <- as.character(merge.SNG$adjustedcelltype)
names(adjustedcelltype) <- names(merge.SNG$adjustedcelltype)

#NK
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), idents = c("NK")),
               feature1 = "cite_CD56", feature2 = "cite_CD16")
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] > 8 & adt["CD16",] < 2.5))
                           ] <- "NK_CD56hiCD16lo"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD16",] > 2.5))
                 ] <- "NK_CD16hi"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] < 8 & adt["CD16",] < 2.5))
                 ] <- "NK_CD56loCD16lo"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort"), idents = c("NK")), 
               # feature1 = "cite_CD56", feature2 = "cite_CD16")
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] > 6 & adt["CD16",] < 1))
                 ] <- "NK_CD56hiCD16lo"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD16",] > 1))
                 ] <- "NK_CD16hi"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] < 6 & adt["CD16",] < 1))
                 ] <- "NK_CD56loCD16lo"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort"), idents = c("NK")), 
#                feature1 = "cite_CD56", feature2 = "cite_CD16")
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] > 5 & adt["CD16",] < 2))
                 ] <- "NK_CD56hiCD16lo"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD16",] > 2))
                 ] <- "NK_CD16hi"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] < 5 & adt["CD16",] < 2))
                 ] <- "NK_CD56loCD16lo"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort"), idents = c("NK")),
               feature1 = "cite_CD56", feature2 = "cite_CD16")
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] > 5 & adt["CD16",] < .5))
                 ] <- "NK_CD56hiCD16lo"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD16",] > .5))
                 ] <- "NK_CD16hi"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] < 5 & adt["CD16",] < .5))
                 ] <- "NK_CD56loCD16lo"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort"), idents = c("NK")), 
#                feature1 = "cite_CD56", feature2 = "cite_CD16")
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] > 8 & adt["CD16",] < .5))
                 ] <- "NK_CD56hiCD16lo"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD16",] > .5))
                 ] <- "NK_CD16hi"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                               merge.SNG$coursecelltype == "NK" & 
                               adt["CD56",] < 8 & adt["CD16",] < .5))
                 ] <- "NK_CD56loCD16lo"
#B cells
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), idents = c("B")),
               feature1 = "cite_IgD", feature2 = "cite_CD27")
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] < 2))
                 ] <- "B_Naive"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDneg"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] < 2))
                 ] <- "B_CD27negIgDneg"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort"), idents = c("B")), 
#                feature1 = "cite_IgD", feature2 = "cite_CD27")
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] < 2))
                 ] <- "B_Naive"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDneg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] < 2))
                 ] <- "B_CD27negIgDneg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort"), idents = c("B")),
               feature1 = "cite_IgD", feature2 = "cite_CD27")
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] > 2.5 & adt["CD27",] < 2))
                 ] <- "B_Naive"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] > 2))
                 ] <- "B_Mem.IgDneg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgD",] < 2.5 & adt["CD27",] < 2))
                 ] <- "B_CD27negIgDneg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort"), idents = c("B")),
               feature1 = "cite_CD11c", feature2 = "cite_IgM")
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] >= 1 & adt["CD11c",] >= 1))
                 ] <- "B_MemSort.IgMpos.CD11cpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] < 1 & adt["CD11c",] < 1))
                 ] <- "B_MemSort.IgMnegCD11cNeg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] < 1 & adt["CD11c",] >= 1))
                 ] <- "B_MemSort.IgMnegCD11cpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] >= 1 & adt["CD11c",] < 1))
                 ] <- "B_MemSort.IgMposCD11cneg"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort"), idents = c("B")), 
#                feature1 = "cite_CD11c", feature2 = "cite_IgM")
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] >= 3 & adt["CD11c",] >= 1))
                 ] <- "B_MemSort.IgMpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] < 3 & adt["CD11c",] < 1))
                 ] <- "B_MemSort.IgMnegCD11cNeg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] < 3 & adt["CD11c",] >= 1))
                 ] <- "B_MemSort.IgMnegCD11cpos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "B" & 
                           adt["IgM",] >= 3 & adt["CD11c",] , 1))
                 ] <- "B_MemSort.IgMposCD11cneg"

#CD4 cells
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), idents = c("CD4")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28")
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] >= 1))
                 ] <- "CD4_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] < 1))
                 ] <- "CD4_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 2 & adt["CD28",] >= 1))
                 ] <- "CD4_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 2 & adt["CD28",] < 1))
                 ] <- "CD4_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort" & cite_CD28 < 10), idents = c("CD4")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1.25 & adt["CD28",] >= .75))
                 ] <- "CD4_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1.25 & adt["CD28",] < .75))
                 ] <- "CD4_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1.25 & adt["CD28",] >= .75))
                 ] <- "CD4_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1.25 & adt["CD28",] < .75))
                 ] <- "CD4_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort" & cite_CD28 < 10), idents = c("CD4")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1 & adt["CD28",] >= .75))
                 ] <- "CD4_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1 & adt["CD28",] < .75))
                 ] <- "CD4_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1 & adt["CD28",] >= .75))
                 ] <- "CD4_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1 & adt["CD28",] < .75))
                 ] <- "CD4_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort" & cite_CD28 < 10), idents = c("CD4")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 0.5 & adt["CD28",] >= 0))
                 ] <- "CD4_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 0.5 & adt["CD28",] < 0))
                 ] <- "CD4_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 0.5 & adt["CD28",] >= 0))
                 ] <- "CD4_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 0.5 & adt["CD28",] < 0))
                 ] <- "CD4_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort" & cite_CD28 < 10), idents = c("CD4")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1 & adt["CD28",] >= 0.75))
                 ] <- "CD4_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] >= 1 & adt["CD28",] < 0.75))
                 ] <- "CD4_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1 & adt["CD28",] >= 0.75))
                 ] <- "CD4_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD4" & 
                           adt["CD45RA",] < 1 & adt["CD28",] < 0.75))
                 ] <- "CD4_CD45RAnegCD28neg"
#CD8 cells
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), idents = c("CD8")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28")
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 3 & adt["CD28",] >= 1))
                 ] <- "CD8_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 3 & adt["CD28",] < 1))
                 ] <- "CD8_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 3 & adt["CD28",] >= 1))
                 ] <- "CD8_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 3 & adt["CD28",] < 1))
                 ] <- "CD8_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort" & cite_CD28 < 10), idents = c("CD8")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] >= .75))
                 ] <- "CD8_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] < .75))
                 ] <- "CD8_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 2 & adt["CD28",] >= .75))
                 ] <- "CD8_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 2 & adt["CD28",] < .75))
                 ] <- "CD8_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort" & cite_CD28 < 10), idents = c("CD8")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 1.5 & adt["CD28",] >= .5))
                 ] <- "CD8_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 1.5 & adt["CD28",] < .5))
                 ] <- "CD8_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 1.5 & adt["CD28",] >= .5))
                 ] <- "CD8_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 1.5 & adt["CD28",] < .5))
                 ] <- "CD8_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort" & cite_CD28 < 10), idents = c("CD8")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 0.75 & adt["CD28",] >= 1))
                 ] <- "CD8_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 0.75 & adt["CD28",] < 1))
                 ] <- "CD8_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 0.75 & adt["CD28",] >= 1))
                 ] <- "CD8_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 0.75 & adt["CD28",] < 1))
                 ] <- "CD8_CD45RAnegCD28neg"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort" & cite_CD28 < 10), idents = c("CD8")),
               feature1 = "cite_CD45RA", feature2 = "cite_CD28", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] >= 0))
                 ] <- "CD8_CD45RAposCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] >= 2 & adt["CD28",] < 0))
                 ] <- "CD8_CD45RAposCD28neg"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 2 & adt["CD28",] >= 0))
                 ] <- "CD8_CD45RAnegCD28pos"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$coursecelltype == "CD8" & 
                           adt["CD45RA",] < 2 & adt["CD28",] < 0))
                 ] <- "CD8_CD45RAnegCD28neg"


#Classical/Intermediate Mono cells
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), subset = mergedcelltype == c("Mono_Class.Int")), 
#                feature1 = "cite_CD14", feature2 = "cite_CD16", pt.size = 0)
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] >= 2))
                 ] <- "Mono_Intermediate"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] < 2))
                 ] <- "Mono_Classical"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort"), subset = mergedcelltype == c("Mono_Class.Int")), 
#                feature1 = "cite_CD14", feature2 = "cite_CD16", pt.size = 0)
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] >= .5))
                 ] <- "Mono_Intermediate"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] < .5))
                 ] <- "Mono_Classical"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort"), subset = mergedcelltype == c("Mono_Class.Int")), 
#                feature1 = "cite_CD14", feature2 = "cite_CD16", pt.size = 0)
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] >= 2.25))
                 ] <- "Mono_Intermediate"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] < 2.25))
                 ] <- "Mono_Classical"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort"), subset = mergedcelltype == c("Mono_Class.Int")), 
#                feature1 = "cite_CD14", feature2 = "cite_CD16", pt.size = 0)
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] >= 0))
                 ] <- "Mono_Intermediate"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] < 0))
                 ] <- "Mono_Classical"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort"), subset = mergedcelltype == c("Mono_Class.Int")), 
#                feature1 = "cite_CD14", feature2 = "cite_CD16", pt.size = 0)
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] >= 0))
                 ] <- "Mono_Intermediate"
adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
                           merge.SNG$mergedcelltype == "Mono_Class.Int" & 
                           adt["CD16",] < 0))
                 ] <- "Mono_Classical"
#pDC/cDC cells
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B1UnSort"), subset = mergedcelltype == c("pDC.cDC")), 
#                feature1 = "cite_CD11c", feature2 = "cite_CD123", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] >= 7.5))
                 ] <- "pDC"
adjustedcelltype[c(which(merge.SNG$Batch == "B1UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] < 7.5))
                 ] <- "cDC"
# FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2UnSort"), subset = mergedcelltype == c("pDC.cDC")), 
#                feature1 = "cite_CD11c", feature2 = "cite_CD123", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] >= 7.5))
                 ] <- "pDC"
adjustedcelltype[c(which(merge.SNG$Batch == "B2UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] < 7.5))
                 ] <- "cDC"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3UnSort"), subset = mergedcelltype == c("pDC.cDC")),
               feature1 = "cite_CD11c", feature2 = "cite_CD123", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] >= 5))
                 ] <- "pDC"
adjustedcelltype[c(which(merge.SNG$Batch == "B3UnSort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] < 5))
                 ] <- "cDC"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B2Sort"), subset = mergedcelltype == c("pDC.cDC")),
               feature1 = "cite_CD11c", feature2 = "cite_CD123", pt.size = 0.1)
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] >= 5))
                 ] <- "pDC"
adjustedcelltype[c(which(merge.SNG$Batch == "B2Sort" & 
                           merge.SNG$mergedcelltype == "pDC.cDC" & 
                           adt["CD123",] < 5))
                 ] <- "cDC"
FeatureScatter(subset(subset(merge.SNG, subset = Batch == "B3Sort"), subset = mergedcelltype == c("pDC.cDC")),
               feature1 = "cite_CD11c", feature2 = "cite_CD123", pt.size = 0.1)
# adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
#                            merge.SNG$mergedcelltype == "pDC.cDC" & 
#                            adt["CD123",] >= 5))
#                  ] <- "pDC"
# adjustedcelltype[c(which(merge.SNG$Batch == "B3Sort" & 
#                            merge.SNG$mergedcelltype == "pDC.cDC" & 
#                            adt["CD123",] < 5))
#                  ] <- "cDC"


write.csv(adjustedcelltype, file = "Metadata/adjustedcelltype.barcodelabelled.20200716.csv")


###############################################
### Re-cluster within the celltypes
###############################################

CITE = as.matrix(GetAssayData(merge.SNG, assay = "CITE", slot = "data"))
limmaBatchNormCITE = limma::removeBatchEffect(CITE, batch = merge.SNG$Batch)
merge.SNG[["limmaCITE"]] <- CreateAssayObject(data = limmaBatchNormCITE)

merge.SNG$batchmajorclusters = paste(as.character(merge.SNG$Batch), as.character(merge.SNG$coursecelltype), sep = "_")
batchmajorclusters = unique(merge.SNG$batchmajorclusters)
# make a list of all the objects with just the major cell types 
batchcelltypeObjects = list()
Idents(merge.SNG) <- "batchmajorclusters"
table(Idents(merge.SNG))
for(i in 1:length(batchmajorclusters)){
  batchcelltypeObjects[[i]] = subset(merge.SNG, idents = batchmajorclusters[i])
}

f1 <- pOverA(0.01, 2)
ffun <- filterfun(f1)
WBADTfilt = list()
# re-cluster within each cell type
for(i in 1:length(batchmajorclusters)){
  WBADTfilt[[i]] = genefilter(GetAssayData(batchcelltypeObjects[[i]][["limmaCITE"]], slot = "data"), ffun)
  adt.dist.temp <- parDist(t(as.matrix(GetAssayData(batchcelltypeObjects[[i]][["limmaCITE"]], slot = "data")))[,which(WBADTfilt[[i]])], threads = 16)
  batchcelltypeObjects[[i]][["rerun_adt_snn"]] <- FindNeighbors(adt.dist.temp, nn.eps = 1)$snn
  batchcelltypeObjects[[i]] <- FindClusters(batchcelltypeObjects[[i]], resolution = c(0.1,0.5,1), graph.name = "rerun_adt_snn", algorithm = 1)
}

clustersperbatchmajorcluster = character()
for(i in 1:length(batchmajorclusters)){
  clustersperbatchmajorcluster[i] = length(unique(batchcelltypeObjects[[i]]$rerun_adt_snn_res.1))
}

clustersperbatchmajorcluster

dir.create("withinBatchWithinCelltypeClustOut_200710", showWarnings = FALSE)
batchaver.obj.1 = list()

### plots of all within-celltype clusters
for(i in 1:length(batchmajorclusters)){
  Idents(batchcelltypeObjects[[i]]) <- "rerun_adt_snn_res.1"
  batchaver.obj.1[[i]] = AverageExpression_MeanOnly(batchcelltypeObjects[[i]], assays = "limmaCITE", return.seurat=F)
  mat_breaks.temp <- quantile_breaks(as.matrix(batchaver.obj.1[[i]]$limmaCITE[WBADTfilt[[i]],]), n = 101)
  ggsave(filename = paste0("withinBatchWithinCelltypeClustOut_200710/",batchmajorclusters[i],"COVID.batchaverExp_withinCelltypeClust_CITE.1.pdf"), width = 8, height = 15,
         plot = pheatmap(batchaver.obj.1[[i]]$limmaCITE[WBADTfilt[[i]],], scale = "none", border_color=NA, viridis(n=100)))
  pheatmap(prop.table(x = table(batchcelltypeObjects[[i]]$rerun_adt_snn_res.1, adjustedcelltype[colnames(batchcelltypeObjects[[i]])]), margin =  2)  %>% '*'(100) %>% round(2),
           filename = paste0("withinBatchWithinCelltypeClustOut_200710/",batchmajorclusters[i],".COVID.withinCelltypeClustPercentageperAdjustedcelltypelabels.1.pdf"), width = 3, height = 6,
           cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
  write.csv(prop.table(x = table(batchcelltypeObjects[[i]]$rerun_adt_snn_res.1, adjustedcelltype[colnames(batchcelltypeObjects[[i]])]), margin =  2)  %>% '*'(100) %>% round(2),
            file = paste0("withinBatchWithinCelltypeClustOut_200710/",batchmajorclusters[i],".COVID.withinCelltypeClustPercentageperAdjustedcelltypelabels.1.csv"))
  pheatmap(prop.table(x = table(batchcelltypeObjects[[i]]$rerun_adt_snn_res.1, paste(batchcelltypeObjects[[i]]$Donor, batchcelltypeObjects[[i]]$Timepoint, sep="_")), margin =  2)  %>% '*'(100) %>% round(2),
           filename = paste0("withinBatchWithinCelltypeClustOut_200710/",batchmajorclusters[i],".COVID.withinCelltypeClustPercentageperDonorTimepoint.1.pdf"), width = 10, height = 10,
           cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = FALSE)
  
  
  ###make a histogram for each cluster at res 1
  # get metadata
  md = batchcelltypeObjects[[i]]@meta.data %>% select(rerun_adt_snn_res.1)
  # get data
  adt = GetAssayData(batchcelltypeObjects[[i]][["limmaCITE"]]) %>%
    t %>%
    as.data.frame %>%
    rownames_to_column("cell")
  md %<>% rownames_to_column("cell")
  md = md[match(x = md$cell, table = adt$cell),  ]
  adt %<>%
    select(-c(cell)) %>%
    mutate(cell_type = md$rerun_adt_snn_res.1) %>%
    select(cell_type, everything())
  
  # get hclust order of proteins and cell types
  x = pheatmap::pheatmap(batchaver.obj.1[[i]]$limmaCITE[WBADTfilt[[i]],],
                         cluster_rows = T, cluster_cols = T,
                         fontsize_col = 10, fontsize_row = 8, border_color = NA)
  dev.off()
  prot_order = rownames(batchaver.obj.1[[i]]$limmaCITE[x$tree_row$order, ])
  celltype_order = colnames(batchaver.obj.1[[i]]$limmaCITE[,x$tree_col$order ])
  
  adt.l = adt %>%
    gather(key = prot, value = normalized_count, CD80:DR3) %>%
    filter(prot %in% prot_order)
  
  adt.l$prot = factor(adt.l$prot)
  adt.l$prot =  reorder.factor(adt.l$prot, new.order = prot_order)
  adt.l$cell_type = factor(as.character(adt.l$cell_type))
  adt.l$cell_type =  reorder.factor(adt.l$cell_type, new.order = celltype_order)
  
  p =  ggplot(adt.l, aes(x=normalized_count, y=prot, fill = ..x..)) +
    geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
    scale_fill_viridis(name = "normalized count", option = "B") +
    scale_x_continuous(limits = c(-1,25)) +
    theme_ridges(font_size = 8, line_size = 0.01) +
    facet_grid(~cell_type)
  
  pdf(paste0("withinBatchWithinCelltypeClustOut_200710/",batchmajorclusters[i],"res.1.protHistogram.pdf"), 24, 16)
  print(p)
  dev.off()
  
}


## Find top protein markers of each cluster

dir.create("DifferentialExpression", showWarnings = FALSE)
clustCITEmarkers = list()
for(i in 1:length(batchmajorclusters)){
  clustCITEmarkers[[i]] = FindAllMarkers(batchcelltypeObjects[[i]], assay = "CITE")
  write.csv(clustCITEmarkers[[i]], file = paste0("DifferentialExpression/", batchmajorclusters[i],"_clusterProteinMarkers.csv"))
}

##relabel BatchCelltypeClusters
for(i in 1:length(batchmajorclusters)){
  batchcelltypeObjects[[i]]$BatchCelltypeClusters <- paste0(batchcelltypeObjects[[i]]$batchmajorclusters, "_", batchcelltypeObjects[[i]]$rerun_adt_snn_res.1)
  Idents(batchcelltypeObjects[[i]]) <- "BatchCelltypeClusters"
  batchaver.obj.1[[i]] = AverageExpression_MeanOnly(batchcelltypeObjects[[i]], assays = "limmaCITE", return.seurat=F)
}

names(batchaver.obj.1) <- batchmajorclusters
remerge.SNG = merge(batchcelltypeObjects[[1]], batchcelltypeObjects[2:length(batchcelltypeObjects)])

#name the clusters. bring back in the clusternames here from csv files
WCTlabels = read.csv(
  file = "withinBatchWithinCelltypeClustOut_200710/withinbatchandcelltypeclustersVsGatingOverlap.labelled3.csv",
  header = TRUE)
Idents(remerge.SNG) <- "BatchCelltypeClusters"
cluster.ids <- WCTlabels$BatchCelltypeClusters
setdiff(Idents(remerge.SNG), cluster.ids)
setdiff(cluster.ids, Idents(remerge.SNG))
remerge.SNG$WCTmergedcelltype <- plyr::mapvalues(x = Idents(remerge.SNG), from = cluster.ids, to = as.character(WCTlabels$WCT.mergedcelltype))
remerge.SNG$WCTcoursecelltype <- plyr::mapvalues(x = Idents(remerge.SNG), from = cluster.ids, to = as.character(WCTlabels$WCT.coursecelltype))

pheatmap(table(remerge.SNG$WCTmergedcelltype, adjustedcelltype[colnames(remerge.SNG)]),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsAdjustedCelltypeOverlap.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", display_numbers = FALSE,
         width = 12, height = 15)

Idents(remerge.SNG) <- "WCTcoursecelltype"

## gating-based cluster QC

QC = rep(TRUE, times = length(colnames(remerge.SNG)))
limmaadt = as.matrix(GetAssayData(remerge.SNG, assay = "limmaCITE", slot = "data"))
FeatureScatter(subset(remerge.SNG, idents = "CD4_Naive"),
               "limmacite_CD4","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Naive" & 
             limmaadt["CD4",] < 1.25))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Naive" &
             limmaadt["CD8",] > 2.75))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Naive" & 
             limmaadt["CD19",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Naive" &
             limmaadt["CD64",] >5))] <- FALSE
FeatureScatter(subset(remerge.SNG, idents = "CD8_Naive"),
               "limmacite_CD4","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Naive" & 
             limmaadt["CD4",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Naive" & 
             limmaadt["CD8",] < 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Naive" & 
             limmaadt["CD19",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Naive" &
             limmaadt["CD64",] > 5))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "CD4_Mem"),
               "limmacite_CD4","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Mem" & 
             limmaadt["CD4",] < 1.25))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Mem" &
             limmaadt["CD8",] > 4))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Mem" & 
             limmaadt["CD19",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD4_Mem" & 
             limmaadt["CD64",] > 5))] <- FALSE
FeatureScatter(subset(remerge.SNG, idents = "CD8_Mem"),
               "limmacite_CD4","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Mem" & 
             limmaadt["CD4",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Mem" &
             limmaadt["CD8",] < 1))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Mem" & 
             limmaadt["CD19",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "CD8_Mem" &
             limmaadt["CD64",] > 5))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "Mono_NonClassical"),
               "limmacite_CD56","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_NonClassical" & 
             limmaadt["CD56",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_NonClassical" &
             limmaadt["CD8",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_NonClassical" & 
             limmaadt["CD19",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_NonClassical" & 
             limmaadt["CD7",] > 2.5))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "Mono_Classical"),
               "limmacite_CD56","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_Classical" & 
             limmaadt["CD56",] > 12))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_Classical" &
             limmaadt["CD8",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_Classical" & 
             limmaadt["CD19",] > 4))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "Mono_Classical" &
             limmaadt["CD7",] > 3))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "B_Naive"),
               "limmacite_CD56","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Naive" & 
             limmaadt["CD56",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Naive" &
             limmaadt["CD8",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Naive" & 
             limmaadt["CD64",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Naive" & 
             limmaadt["CD4",] > 2))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "B_Mem"),
               "limmacite_CD56","limmacite_CD8")
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Mem" & 
             limmaadt["CD56",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Mem" & 
             limmaadt["CD8",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Mem" & 
             limmaadt["CD64",] > 5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "B_Mem" &
             limmaadt["CD4",] > 2.5))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "NK_CD16hi"),
               "limmacite_CD56","limmacite_CD19")
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD16hi" & 
             limmaadt["CD19",] > 3))
   ] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD16hi" & 
             limmaadt["CD64",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD16hi" &
             limmaadt["CD4",] > 2))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "NK_CD56loCD16lo"),
               "limmacite_CD56","limmacite_CD19")
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56loCD16lo" &
             limmaadt["CD19",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56loCD16lo" & 
             limmaadt["CD64",] > 3 ))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56loCD16lo" &
             limmaadt["CD4",] > 2))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "NK_CD56hiCD16lo"),
               "limmacite_CD56","limmacite_CD19")
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56hiCD16lo" &
             limmaadt["CD19",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56hiCD16lo" & 
             limmaadt["CD64",] > 3 ))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "NK_CD56hiCD16lo" &
             limmaadt["CD4",] > 2))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "pDC"),
               "limmacite_CD123","limmacite_CD19")
QC[c(which(remerge.SNG$WCTcoursecelltype == "pDC" & 
             limmaadt["CD123",] < 4 ))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "pDC" &
             limmaadt["CD19",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "pDC" &
             limmaadt["CD16",] > 2.5))] <- FALSE


FeatureScatter(subset(remerge.SNG, idents = "cDC"),
               "limmacite_CD11c","limmacite_CD19")
QC[c(which(remerge.SNG$WCTcoursecelltype == "cDC" & 
             limmaadt["CD8",] > 2 ))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "cDC" &
             limmaadt["CD19",] > 3))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "cDC" &
             limmaadt["CD7",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "cDC" &
             limmaadt["CD16",] > 2))] <- FALSE

FeatureScatter(subset(remerge.SNG, idents = "PB_Plasmablasts"),
               "limmacite_CD38","limmacite_CD23")
QC[c(which(remerge.SNG$WCTcoursecelltype == "PB_Plasmablasts" & 
             limmaadt["CD8",] > 2 ))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "PB_Plasmablasts" &
             limmaadt["CD4",] > 2))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "PB_Plasmablasts" &
             limmaadt["CD11b",] > 2.5))] <- FALSE
QC[c(which(remerge.SNG$WCTcoursecelltype == "PB_Plasmablasts" &
             limmaadt["CD16",] > 2.5))] <- FALSE

WCT.Labels = data.frame(WCTcoursecelltype = remerge.SNG$WCTcoursecelltype, 
                        WCTmergedcelltype = remerge.SNG$WCTmergedcelltype, 
                        QC = QC,
                        coursecelltype = remerge.SNG$coursecelltype, 
                        adjustedcelltype = adjustedcelltype[colnames(remerge.SNG)])
write.csv(WCT.Labels, file = "withinBatchWithinCelltypeClustOut_200710/WCT.Labels.csv")

pheatmap(table(remerge.SNG$WCTcoursecelltype[QC], adjustedcelltype[colnames(remerge.SNG)[QC]]),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsAdjustedCelltypeOverlap.Course_QCd.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", display_numbers = FALSE,
         width = 12, height = 15)
pheatmap(table(remerge.SNG$WCTmergedcelltype[QC], adjustedcelltype[colnames(remerge.SNG)[QC]]),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsAdjustedCelltypeOverlap_QCd.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", display_numbers = FALSE,
         width = 12, height = 15)
pheatmap(prop.table(x = table(remerge.SNG$WCTmergedcelltype[QC], paste(remerge.SNG$Donor[QC], remerge.SNG$Timepoint[QC], sep="_")), margin = 2),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsDonor_Timepoint.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "none", display_numbers = FALSE,
         width = 12, height = 15)
pheatmap(prop.table(x = table(remerge.SNG$WCTmergedcelltype[QC], paste(remerge.SNG$Donor[QC], remerge.SNG$Timepoint[QC], sep="_")), margin = 2),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsDonor_Timepoint.rowscaled.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", display_numbers = FALSE,
         width = 12, height = 15)
pheatmap(log1p(prop.table(x = table(remerge.SNG$WCTmergedcelltype[QC], paste(remerge.SNG$Donor[QC], remerge.SNG$Timepoint[QC], sep="_")), margin = 2)),
         filename = "withinBatchWithinCelltypeClustOut_200710/WCTclustVsDonor_Timepoint.log1p.png",
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "none", display_numbers = FALSE,
         width = 12, height = 15)

write.csv(table(remerge.SNG$WCTmergedcelltype[QC], paste(remerge.SNG$Donor[QC], remerge.SNG$Timepoint[QC], sep="_")), 
          file = "withinBatchWithinCelltypeClustOut_200710/WCTclustCellsPerDonor_Timepoint.csv")
write.csv(table(remerge.SNG$WCTcoursecelltype[QC], paste(remerge.SNG$Donor[QC], remerge.SNG$Timepoint[QC], sep="_")), 
          file = "withinBatchWithinCelltypeClustOut_200710/WCTclustCellsPerDonor_Timepoint.Course.csv")
write.csv(table(remerge.SNG$WCTcoursecelltype[QC], remerge.SNG$Batch[QC]), 
          file = "withinBatchWithinCelltypeClustOut_200710/WCTclustCellsPerBatch.Course.csv")
write.csv(table(remerge.SNG$WCTmergedcelltype[QC], remerge.SNG$Batch[QC]), 
          file = "withinBatchWithinCelltypeClustOut_200710/WCTclustCellsPerBatch.csv")

saveRDS(remerge.SNG, "SeuratObjects/remerge.SNG.withWCTclusters.Rds")

