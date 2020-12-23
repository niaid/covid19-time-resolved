### normalize and parse pseudobulk for plotting
### the input is the output data from pseudobulk object generation pipeline
library(Biobase)
library(edgeR)
library(plyr)
library(tidyverse)
### unsorted #####################################################################################
DGE_IN_PATH <- "input/differential_expression/"
# UnSorted.WCTcourse <- readRDS(paste(DGE_IN_PATH, "all_samples/pseudobulk_dgelists_unfiltered/Unsorted-WCTcoursecelltype.rds", sep = ""))
UnSorted.WCTcourse <- readRDS(paste0(DGE_IN_PATH, "/Unsorted-WCTcoursecelltype.rds"))

eset_list <- lapply(UnSorted.WCTcourse, function(dge){
  mat <- cpm(dge$counts, log = TRUE)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})

# parse phenodata
# also add steroid data
steroid.meta <- read.csv("input/batch.citeseq.samples.steroid.csv", header = TRUE, stringsAsFactors = FALSE)
steroid.meta$sample_id <- paste(steroid.meta$subject_id, steroid.meta$visit, sep = "_")

parse_meta <- function(dge){
  dge@phenoData@data$PC1_cat <- as.character(dge@phenoData@data$PC1_cat)
  dge@phenoData@data$PC1_cat <- replace(dge@phenoData@data$PC1_cat, 
                                        dge@phenoData@data$Class == "HC", "HC")
  dge@phenoData@data$PC1_cat <- factor(dge@phenoData@data$PC1_cat, levels = c("HC","PC1_low","PC1_high"))
  dge$Subject <- sapply(str_split(rownames(dge@phenoData@data), pattern = "_"), function(x)x[1])
  dge$sample_id <- paste(dge$Subject, dge$Timepoint, sep = "_")
  dge$steroid.use <- plyr::mapvalues(x = dge$sample_id, from = steroid.meta$sample_id, to = steroid.meta$steroid.use)
  dge@phenoData@data$steroid.use <- replace(dge@phenoData@data$steroid.use, 
                                            dge@phenoData@data$Class == "HC", "HC")
  dge$steroid.use2 <- plyr::mapvalues(x = dge$sample_id, from = steroid.meta$sample_id, to = steroid.meta$steroid.use2)
  dge@phenoData@data$steroid.use2 <- replace(dge@phenoData@data$steroid.use2, 
                                             dge@phenoData@data$Class == "HC", "HC")
  dge$o2.supplement <- plyr::mapvalues(x = dge$sample_id, from = steroid.meta$sample_id, to = steroid.meta$o2.supplement)
  dge@phenoData@data$o2.supplement <- replace(dge@phenoData@data$o2.supplement, 
                                              dge@phenoData@data$Class == "HC", "HC")
  # filter <- colnames(dge)[!is.na(dge@phenoData@data$PC1_cat)]
  dge <- dge
}

eset_list <- lapply(eset_list, parse_meta)

dir.create("output/dge_lists")
saveRDS(eset_list, "output/dge_lists/pbulk_eset_list_normalized_WCTcourse_metafiltered.rds")

### sorted #####################################################################################
# Sorted.specificgate <- readRDS(paste(DGE_IN_PATH, "all_samples/pseudobulk_dgelists_unfiltered/Sorted-specific_gating.rds", sep = ""))
Sorted.specificgate <- readRDS(paste0(DGE_IN_PATH, "/Sorted-specific_gating.rds"))

eset_list_sorted <- lapply(Sorted.specificgate, function(dge){
  mat <- cpm(dge$counts, log = TRUE)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})

# parse phenodata
# also add steroid data
eset_list_sorted <- lapply(eset_list_sorted, parse_meta)

saveRDS(eset_list_sorted, "output/dge_lists/pbulk_eset_list_normalized_specific_gating_metafiltered.rds")


