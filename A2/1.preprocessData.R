#this is tested using R 3.6.1 on a high-performance comupting node with 16 cores and at least 160 gb of ram. 
library("Seurat") #load Seurat 3.2.2
library("dplyr")
library("matrixStats")
library('tidyverse')

#define 10x lanes for import
B1Lanes = c(1:4,5,7,8)
B2LanesUnsorted = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16")
B2LanesSorted = c("17","18","19","20","21","22","23","24")
B3LanesUnsorted = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16")
B3LanesSorted = c("17","18","19","20","21","22","23","24")

#initialize objects for the data
B1_data = list()
B1_SeuratObj = list()
B2_US_data = list()
B2_S_data = list()
B2_US_SeuratObj = list()
B2_S_SeuratObj = list()
B3_US_data = list()
B3_S_data = list()
B3_US_SeuratObj = list()
B3_S_SeuratObj = list()

### Import the data
# the code below assumes the HDF5 data has been downloaded from GEO into a directory labelled "Counts"
for(i in 1:length(B1Lanes)){
  B1_data[[i]] = list(RNA = Read10X_h5(paste("Counts/B1_10xlane",B1Lanes[[i]],"_RNA_filtered_feature_bc_matrix.h5", sep="")),
                      CSP = Read10X_h5(paste("Counts/B1_10xlane",B1Lanes[[i]],"_CSP_filtered_feature_bc_matrix.h5", sep=""))
  )
  sharedBC = intersect(colnames(B1_data[[i]]$RNA$`Gene Expression`), colnames(B1_data[[i]]$CSP))
  B1_SeuratObj[[i]] <- CreateSeuratObject(counts = B1_data[[i]]$RNA$'Gene Expression'[,sharedBC], assay = "RNA", min.feature = 5)
  B1_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_data[[i]]$CSP[11:202,sharedBC])
  B1_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B1_data[[i]]$CSP[c(1:3,5),sharedBC])
  B1_SeuratObj[[i]] <- RenameCells(B1_SeuratObj[[i]], new.names = paste(substr(colnames(B1_SeuratObj[[i]]), start = 1, stop = 17),B1Lanes[[i]], sep = ""))
  B1_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_SeuratObj[[i]])))
}
for(i in 1:length(B2LanesUnsorted)){
  B2_US_data[[i]] = list(RNA = Read10X_h5(paste("Counts/B2_10xlane",B2LanesUnsorted[[i]],"_RNA_filtered_feature_bc_matrix.h5", sep="")),
                         CSP = Read10X_h5(paste("Counts/B2_10xlane",B2LanesUnsorted[[i]],"_CSP_filtered_feature_bc_matrix.h5", sep=""))
  )
  sharedBC = intersect(colnames(B2_US_data[[i]]$RNA$`Gene Expression`), colnames(B2_US_data[[i]]$CSP))
  B2_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B2_US_data[[i]]$RNA$'Gene Expression'[,sharedBC], assay = "RNA", min.feature = 5)
  B2_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B2_US_data[[i]]$CSP[11:202,sharedBC])
  B2_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B2_US_data[[i]]$CSP[c(1:3,5),sharedBC])
  B2_US_SeuratObj[[i]] <- RenameCells(B2_US_SeuratObj[[i]], new.names = paste(substr(colnames(B2_US_SeuratObj[[i]]), start = 1, stop = 17),B2LanesUnsorted[[i]], sep = ""))
  B2_US_SeuratObj[[i]]$Batch  <- rep("B2UnSort", length(colnames(B2_US_SeuratObj[[i]])))
}
for(i in 1:length(B2LanesSorted)){
  B2_S_data[[i]] = list(RNA = Read10X_h5(paste("Counts/B2_10xlane",B2LanesSorted[[i]],"_RNA_filtered_feature_bc_matrix.h5", sep="")),
                        CSP = Read10X_h5(paste("Counts/B2_10xlane",B2LanesSorted[[i]],"_CSP_filtered_feature_bc_matrix.h5", sep=""))
  )
  sharedBC = intersect(colnames(B2_S_data[[i]]$RNA$`Gene Expression`), colnames(B2_S_data[[i]]$CSP))
  B2_S_SeuratObj[[i]] <- CreateSeuratObject(counts = B2_S_data[[i]]$RNA$'Gene Expression'[,sharedBC], assay = "RNA", min.feature = 5)
  B2_S_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B2_S_data[[i]]$CSP[11:202,sharedBC])
  B2_S_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B2_S_data[[i]]$CSP[c(1:3,5),sharedBC])
  B2_S_SeuratObj[[i]] <- RenameCells(B2_S_SeuratObj[[i]], new.names = paste(substr(colnames(B2_S_SeuratObj[[i]]), start = 1, stop = 17),B2LanesSorted[[i]], sep = ""))
  B2_S_SeuratObj[[i]]$Batch  <- rep("B2Sort", length(colnames(B2_S_SeuratObj[[i]])))
}


for(i in 1:length(B3LanesUnsorted)){
  B3_US_data[[i]] = list(RNA = Read10X_h5(paste("Counts/B3_10xlane",B3LanesUnsorted[[i]],"_RNA_filtered_feature_bc_matrix.h5", sep="")),
                         CSP = Read10X_h5(paste("Counts/B3_10xlane",B3LanesUnsorted[[i]],"_CSP_filtered_feature_bc_matrix.h5", sep=""))
  )
  sharedBC = intersect(colnames(B3_US_data[[i]]$RNA$`Gene Expression`), colnames(B3_US_data[[i]]$CSP))
  B3_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B3_US_data[[i]]$RNA$'Gene Expression'[,sharedBC], assay = "RNA", min.feature = 5)
  B3_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B3_US_data[[i]]$CSP[11:202,sharedBC])
  B3_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B3_US_data[[i]]$CSP[c(1:3,5),sharedBC])
  B3_US_SeuratObj[[i]] <- RenameCells(B3_US_SeuratObj[[i]], new.names = paste(substr(colnames(B3_US_SeuratObj[[i]]), start = 1, stop = 17),B3LanesUnsorted[[i]], sep = ""))
  B3_US_SeuratObj[[i]]$Batch  <- rep("B3UnSort", length(colnames(B3_US_SeuratObj[[i]])))
}
for(i in 1:length(B3LanesSorted)){
  B3_S_data[[i]] = list(RNA = Read10X_h5(paste("Counts/B3_10xlane",B3LanesSorted[[i]],"_RNA_filtered_feature_bc_matrix.h5", sep="")),
                        CSP = Read10X_h5(paste("Counts/B3_10xlane",B3LanesSorted[[i]],"_CSP_filtered_feature_bc_matrix.h5", sep=""))
  )
  sharedBC = intersect(colnames(B3_S_data[[i]]$RNA$`Gene Expression`), colnames(B3_S_data[[i]]$CSP))
  B3_S_SeuratObj[[i]] <- CreateSeuratObject(counts = B3_S_data[[i]]$RNA$'Gene Expression'[,sharedBC], assay = "RNA", min.feature = 5)
  B3_S_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B3_S_data[[i]]$CSP[11:202,sharedBC])
  B3_S_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B3_S_data[[i]]$CSP[c(1:3,5),sharedBC])
  B3_S_SeuratObj[[i]] <- RenameCells(B3_S_SeuratObj[[i]], new.names = paste(substr(colnames(B3_S_SeuratObj[[i]]), start = 1, stop = 17),B3LanesSorted[[i]], sep = ""))
  B3_S_SeuratObj[[i]]$Batch  <- rep("B3Sort", length(colnames(B3_S_SeuratObj[[i]])))
}

names(B1_data) = names(B1_SeuratObj) = B1Lanes
names(B2_US_data) = names(B2_US_SeuratObj) = B2LanesUnsorted
names(B3_US_data) = names(B3_US_SeuratObj) = B3LanesUnsorted
names(B2_S_data) = names(B2_S_SeuratObj) = B2LanesSorted
names(B3_S_data) = names(B3_S_SeuratObj) = B3LanesSorted

# use HTOdemux to identify negatives droplets to use for DSB norm
# merge objects
B1_merge = merge(B1_SeuratObj[[1]], B1_SeuratObj[2:length(B1Lanes)])
B2_US_merge = merge(B2_US_SeuratObj[[1]], B2_US_SeuratObj[2:length(B2LanesUnsorted)])
B3_US_merge = merge(B3_US_SeuratObj[[1]], B3_US_SeuratObj[2:length(B3LanesUnsorted)])
B2_S_merge = merge(B2_S_SeuratObj[[1]], B2_S_SeuratObj[2:length(B2LanesSorted)])
B3_S_merge = merge(B3_S_SeuratObj[[1]], B3_S_SeuratObj[2:length(B3LanesSorted)])

# remove abmormally high count cells
B1_merge <- subset(B1_merge, subset = nCount_CITE < 30000)
B2_US_merge <- subset(B2_US_merge, subset = nCount_CITE < 30000)
B3_US_merge <- subset(B3_US_merge, subset = nCount_CITE < 30000)
B2_S_merge <- subset(B2_S_merge, subset = nCount_CITE < 30000)
B3_S_merge <- subset(B3_S_merge, subset = nCount_CITE < 30000)

B1_merge <- subset(B1_merge, subset = nCount_HTO < 15001)
B2_US_merge <- subset(B2_US_merge, subset = nCount_HTO < 15001)
B3_US_merge <- subset(B3_US_merge, subset = nCount_HTO < 15001)
B2_S_merge <- subset(B2_S_merge, subset = nCount_HTO < 15001)
B3_S_merge <- subset(B3_S_merge, subset = nCount_HTO < 15001)

B1_merge <- NormalizeData(B1_merge, assay = "HTO", normalization.method = "CLR")
B1_merge <- ScaleData(B1_merge, assay = "HTO", model.use = "linear")
B1_merge = HTODemux(B1_merge, positive.quantile = 0.99)
Idents(B1_merge) <- "hash.ID"

B2_US_merge <- NormalizeData(B2_US_merge, assay = "HTO", normalization.method = "CLR")
B2_US_merge <- ScaleData(B2_US_merge, assay = "HTO", model.use = "linear")
B2_US_merge = HTODemux(B2_US_merge, positive.quantile = 0.99)
Idents(B2_US_merge) <- "hash.ID"

B2_S_merge <- NormalizeData(B2_S_merge, assay = "HTO", normalization.method = "CLR")
B2_S_merge <- ScaleData(B2_S_merge, assay = "HTO", model.use = "linear")
B2_S_merge = HTODemux(B2_S_merge, positive.quantile = 0.99)
Idents(B2_S_merge) <- "hash.ID"

B3_US_merge <- NormalizeData(B3_US_merge, assay = "HTO", normalization.method = "CLR")
B3_US_merge <- ScaleData(B3_US_merge, assay = "HTO", model.use = "linear")
B3_US_merge = HTODemux(B3_US_merge, positive.quantile = 0.99)
Idents(B3_US_merge) <- "hash.ID"

B3_S_merge <- NormalizeData(B3_S_merge, assay = "HTO", normalization.method = "CLR")
B3_S_merge <- ScaleData(B3_S_merge, assay = "HTO", model.use = "linear")
B3_S_merge = HTODemux(B3_S_merge, positive.quantile = 0.99)
Idents(B3_S_merge) <- "hash.ID"

#bring in demuxlet SNP-based demultiplexing calls

B1_demuxbestList = list()
for(i in 1:length(B1Lanes)){
  B1_demuxbestList[[i]] = read.table(paste("demuxletResults/B1_10xlane",B1Lanes[[i]],"_demuxed.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length(B1Lanes)){
  B1_demuxbestList[[i]]$NewBarcode = paste(substr(B1_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B1Lanes[[i]], sep = "")
}
B1_demuxbestdf <- plyr::ldply(B1_demuxbestList, data.frame)
length(which(colnames(B1_merge) %in% B1_demuxbestdf$NewBarcode))
setdiff(colnames(B1_merge), B1_demuxbestdf$NewBarcode)
rownames(B1_demuxbestdf) <- B1_demuxbestdf$NewBarcode
B1_merge <- AddMetaData(B1_merge, metadata = B1_demuxbestdf[colnames(B1_merge),])

B2_US_demuxbestList = list()
for(i in 1:length(B2LanesUnsorted)){
  B2_US_demuxbestList[[i]] = read.table(paste("demuxletResults/B2_10xlane",B2LanesUnsorted[[i]],"_demuxed.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length(B2LanesUnsorted)){
  B2_US_demuxbestList[[i]]$NewBarcode = paste(substr(B2_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B2LanesUnsorted[[i]], sep = "")
}
B2_US_demuxbestdf <- plyr::ldply(B2_US_demuxbestList, data.frame)
length(which(colnames(B2_US_merge) %in% B2_US_demuxbestdf$NewBarcode))
setdiff(colnames(B2_US_merge), B2_US_demuxbestdf$NewBarcode)
rownames(B2_US_demuxbestdf) <- B2_US_demuxbestdf$NewBarcode
B2_US_merge <- AddMetaData(B2_US_merge, metadata = B2_US_demuxbestdf[colnames(B2_US_merge),])

B2_S_demuxbestList = list()
for(i in 1:length(B2LanesSorted)){
  B2_S_demuxbestList[[i]] = read.table(paste("demuxletResults/B2_10xlane",B2LanesSorted[[i]],"_demuxed.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length(B2LanesSorted)){
  B2_S_demuxbestList[[i]]$NewBarcode = paste(substr(B2_S_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B2LanesSorted[[i]], sep = "")
}
B2_S_demuxbestdf <- plyr::ldply(B2_S_demuxbestList, data.frame)
length(which(colnames(B2_S_merge) %in% B2_S_demuxbestdf$NewBarcode))
length(setdiff(colnames(B2_S_merge), B2_S_demuxbestdf$NewBarcode))
rownames(B2_S_demuxbestdf) <- B2_S_demuxbestdf$NewBarcode
B2_S_merge <- AddMetaData(B2_S_merge, metadata = B2_S_demuxbestdf[colnames(B2_S_merge),])

B3_US_demuxbestList = list()
for(i in 1:length(B3LanesUnsorted)){
  B3_US_demuxbestList[[i]] = read.table(paste("demuxletResults/B3_10xlane",B3LanesUnsorted[[i]],"_demuxed.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length(B3LanesUnsorted)){
  B3_US_demuxbestList[[i]]$NewBarcode = paste(substr(B3_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B3LanesUnsorted[[i]], sep = "")
}
B3_US_demuxbestdf <- plyr::ldply(B3_US_demuxbestList, data.frame)
length(which(colnames(B3_US_merge) %in% B3_US_demuxbestdf$NewBarcode))
setdiff(colnames(B3_US_merge), B3_US_demuxbestdf$NewBarcode)
rownames(B3_US_demuxbestdf) <- B3_US_demuxbestdf$NewBarcode
B3_US_merge <- AddMetaData(B3_US_merge, metadata = B3_US_demuxbestdf[colnames(B3_US_merge),])

B3_S_demuxbestList = list()
for(i in 1:length(B3LanesSorted)){
  B3_S_demuxbestList[[i]] = read.table(paste("demuxletResults/B3_10xlane",B3LanesSorted[[i]],"_demuxed.best", sep = ""), sep = "\t", header = TRUE)
}
for(i in 1:length(B3LanesSorted)){
  B3_S_demuxbestList[[i]]$NewBarcode = paste(substr(B3_S_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B3LanesSorted[[i]], sep = "")
}
B3_S_demuxbestdf <- plyr::ldply(B3_S_demuxbestList, data.frame)
length(which(colnames(B3_S_merge) %in% B3_S_demuxbestdf$NewBarcode))
length(setdiff(colnames(B3_US_merge), B3_S_demuxbestdf$NewBarcode))
rownames(B3_S_demuxbestdf) <- B3_S_demuxbestdf$NewBarcode
B3_S_merge <- AddMetaData(B3_S_merge, metadata = B3_S_demuxbestdf[colnames(B3_S_merge),])

# merge the data into one object
merge.SNG <- merge(subset(B1_merge, subset = DROPLET.TYPE == "SNG"),list(
  subset(B2_US_merge, subset = DROPLET.TYPE == "SNG"),
  subset(B2_S_merge, subset = DROPLET.TYPE == "SNG"),
  subset(B3_US_merge, subset = DROPLET.TYPE == "SNG"),
  subset(B3_S_merge, subset = DROPLET.TYPE == "SNG")
))
merge.NEG <- merge(subset(B1_merge, subset = DROPLET.TYPE == "AMB" & hash.ID == "Negative"),list(
  subset(B2_US_merge, subset = DROPLET.TYPE == "AMB" & hash.ID == "Negative"),
  subset(B2_S_merge, subset = DROPLET.TYPE == "AMB" & hash.ID == "Negative"),
  subset(B3_US_merge, subset = DROPLET.TYPE == "AMB" & hash.ID == "Negative"),
  subset(B3_S_merge, subset = DROPLET.TYPE == "AMB" & hash.ID == "Negative")
))

#filter low/high umi cells
mito.genes = grep(pattern = "^MT-", x = rownames(merge.NEG), value = TRUE)
(mean(merge.NEG$nFeature_RNA) + 3*sd(merge.NEG$nFeature_RNA)) 
percent.mito.merge.SNG = Matrix::colSums(merge.SNG[mito.genes,])/Matrix::colSums(merge.SNG)
merge.SNG <- AddMetaData(object = merge.SNG, metadata = percent.mito.merge.SNG, col.name = "percent.mito")
VlnPlot(merge.SNG, features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)
merge.SNG = subset(merge.SNG, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mito < 0.30 & nCount_CITE < 15000 & nCount_HTO < 5000)
VlnPlot(merge.SNG, features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)


# # DSB normalization on the HTO
isotype.control.name.vec = c("IgG1Kiso", "IgG2aKiso", "IgG2bKiso", "ratIgG2bKiso" )
source("utilityFunctions/dsb_normalization_functions.R")

norm.HTO.list.merge.SNG = DSBNormalizeProtein(cell.columns.protein.matrix = as.matrix(GetAssayData(merge.SNG[["HTO"]], slot = "counts")),
                                                    control.protein.matrix = as.matrix(GetAssayData(merge.NEG[["HTO"]], slot = "counts")),
                                                    define.pseudocount = TRUE, pseudocount.use = 10, denoise_counts = FALSE,
                                                    isotype.control.name.vec =  NULL)
merge.SNG[["HTO"]] <- SetAssayData(merge.SNG[["HTO"]], slot = "data", new.data = norm.HTO.list.merge.SNG)

 
# ##### manual thresholding hash calls, Unsorted cells
HTONames = rownames(merge.SNG[["HTO"]])

merge.SNG <- ScaleData(merge.SNG, assay = "HTO", vars.to.regress = "Batch")

## Export QC plots to use for manual thresholding 
pdf("QCplots/HTOdist.check.merge.SNG.pdf")
for (i in 1:4) {
  hist(as.matrix(GetAssayData(merge.SNG[["HTO"]], slot = "data"))[i,], breaks=100)
}
dev.off()

# manual threshold HTO calls
merge.SNG.hithres = c(2.25, 3, 2.25, 2)
merge.SNG.lothres = c(2, 2.75, 2, 1.75)
merge.SNG.automaticHTOsHiThres = character(length=length(colnames(merge.SNG)))
for(i in 1:length(colnames(merge.SNG))){
  merge.SNG.automaticHTOsHiThres[i] = paste(HTONames[GetAssayData(merge.SNG[["HTO"]], slot = "data")[,i] > merge.SNG.hithres], collapse="+")
}

# manual lower threshold HTO calls, for doublet removal
merge.SNG.automaticHTOsloThres = character(length=length(colnames(merge.SNG)))
for(i in 1:length(colnames(merge.SNG))){
  merge.SNG.automaticHTOsloThres[i] = paste(HTONames[GetAssayData(merge.SNG[["HTO"]], slot = "data")[,i] > merge.SNG.lothres], collapse="+")
}

merge.SNG.singlets = merge.SNG.automaticHTOsHiThres %in% HTONames & merge.SNG.automaticHTOsloThres %in% HTONames
merge.SNG.autoHashcalls = ifelse(merge.SNG.singlets, merge.SNG.automaticHTOsHiThres, "NonSinglet")
merge.SNG <- AddMetaData(object = merge.SNG, metadata = merge.SNG.autoHashcalls, col.name = "autoHashcalls")

Idents(merge.SNG) <- "autoHashcalls"
FeatureScatter(subset(merge.SNG, idents = c("HTO1","HTO2", "HTO3", "HTO5")), feature1="hto_HTO1", feature2="hto_HTO2")
FeatureScatter(subset(merge.SNG, idents = c("HTO1","HTO2", "HTO3", "HTO5")), feature1="hto_HTO3", feature2="hto_HTO5")
table(merge.SNG$hash.ID, merge.SNG$autoHashcalls)
table(merge.SNG$BEST.GUESS, merge.SNG$autoHashcalls)
table(merge.SNG$BEST.GUESS, merge.SNG$hash.ID)
length(which(!merge.SNG$autoHashcalls == "NonSinglet"))
merge.NonSinglet <- subset(merge.SNG, idents = "NonSinglet")
merge.SNG <- subset(merge.SNG, idents = "NonSinglet", invert = TRUE)

# DSB Normalization
limma.norm.ADT.list.merge.SNG = DSBNormalizeProtein(cell.columns.protein.matrix = as.matrix(GetAssayData(merge.SNG[["CITE"]], slot = "counts")),
                                                    control.protein.matrix = as.matrix(GetAssayData(merge.NEG[["CITE"]], slot = "counts")),
                                                    define.pseudocount = TRUE, pseudocount.use = 10, denoise_counts = TRUE,
                                                    isotype.control.name.vec =  isotype.control.name.vec)



merge.SNG[["CITE"]] <- SetAssayData(merge.SNG[["CITE"]], slot = "data", new.data = limma.norm.ADT.list.merge.SNG$denoised_adt)
merge.SNG <- AddMetaData(object = merge.SNG, metadata = limma.norm.ADT.list.merge.SNG$cellwise_background_mean, col.name = "ADTmclust1mean")


ggsave("QCplots/vlnplot.selectedmarkers.byBatch.merge.SNG.pdf", width = 12, height = 12,
       VlnPlot(merge.SNG, features = c("cite_CD4", "cite_CD3","cite_CD36","cite_CD5","cite_KLRG1","cite_CD278"), 
               group.by = "Batch", pt.size = 0, slot = "data")
)
ggsave("QCplots/vlnplot.selectedmarkers2.byBatch.merge.SNG.pdf", width = 12, height = 12,
       VlnPlot(merge.SNG, features = c("cite_HLA-ABC", "cite_CD11c","cite_CD33","cite_CD19","cite_HLA-DR","cite_CD11b"), 
               group.by = "Batch", pt.size = 0, slot = "data")
)

## add metadata
merge.SNG$Donor = sapply(strsplit(as.character(merge.SNG$BEST.GUESS),split = ","),'[',1)
merge.SNG$Donor = sapply(strsplit(as.character(merge.SNG$Donor),split = "_"),'[',1)
merge.SNG$Sample = paste(merge.SNG$Batch, merge.SNG$Donor,
                         merge.SNG$autoHashcalls, sep = "_")

Idents(merge.SNG) <- "Sample"
table(Idents(merge.SNG))

#save objects for clustering
saveRDS(merge.SNG, file = "SeuratObjects/merge.SNG.forClust.rds")


