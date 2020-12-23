library(edgeR)
library(GSVA)
library(Biobase)
library(reshape2)
library(BiocParallel)
library(readr)
library(tidyverse)
source("util_fun/external.data.gsvadftoelist.R")

### get GSVA score of the filtered pseudobulk objects ##############################################
### correlation with severity -- *covid samples only*
### using LE genes of severe-mild, Schulte-Schrepping et al, 2020, Cell
### input is output from filtered pseudobulk objects
### generate GSVA score files for NK.modulescore.corr.Rmd, output files also provided in input/SchulteSchrepping/cohort2
### cohort2 #######################################################################################
sample.groups <- "all_timepoints_covid_only"
DGELISTS_IN_PATH <- "input/SchulteSchrepping/cohort2/pseudobulk_dgelists_normalized/"

# get filtered data(cell #/sample and genes filtered)
files <- list.files(DGELISTS_IN_PATH, full.names = TRUE)
OUT_DIR <- "input/SchulteSchrepping/cohort2"

genesets <- readRDS("input/kegg_go_btm_reactome_foointerferon.rds")

pbulk_list_cohort2 <- lapply(files, readRDS)
names(pbulk_list_cohort2) <- gsub("-pseudobulk_dgelist_normalized\\.rds", "", basename(files))

# list of expression data by cell type
eset_list_cohort2 <- lapply(pbulk_list_cohort2, function(dge){
  mat <- DESeq2::varianceStabilizingTransformation(dge$counts)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})
eset_list_cohort2 <- unlist(eset_list_cohort2)


# get fgsea results of LE genes from the severity-mild model
FGSEA_IN_PATH <- file.path("input/SchulteSchrepping/cohort2/fgsea_tables/severe-mild")
fgsea_files_cohort2 <- list.files(FGSEA_IN_PATH, full.names = TRUE)

fgsea_list_cohort2 <- lapply(fgsea_files_cohort2, function(path){
  read_tsv(path,
           col_types = cols(
             pathway = col_character(),
             pval = col_double(),
             padj = col_double(),
             ES = col_double(),
             NES = col_double(),
             nMoreExtreme = col_double(),
             size = col_double(),
             leadingEdge = col_character()
           ))
})

names(fgsea_list_cohort2) <- gsub("--model.+\\.tsv", "", basename(fgsea_files_cohort2))
fgsea_res_cohort2 <- bind_rows(fgsea_list_cohort2, .id = "celltype")

scores_list_cohort2 <- lapply(1:length(eset_list_cohort2), function(eset){
  cat(names(eset_list_cohort2)[eset],"\n")
  celltype.gsea.res <- subset(fgsea_res_cohort2, celltype == names(eset_list_cohort2)[eset])
  if (nrow(celltype.gsea.res) > 0) {
    celltype.leading.edge.genesets.df <- aggregate(leadingEdge ~ pathway,celltype.gsea.res,paste,collapse=" ")
    celltype.leading.edge.genesets <- sapply(celltype.leading.edge.genesets.df$leadingEdge,function(x){unique(unlist(strsplit(x," ")))})
    names(celltype.leading.edge.genesets) <- celltype.leading.edge.genesets.df$pathway
    module.scores <- gsva(expr = eset_list_cohort2[[eset]], gset.idx.list = celltype.leading.edge.genesets, method = "gsva", parallel.sz = 16)
    return(cbind(reshape2::melt(exprs(module.scores)),celltype= names(eset_list_cohort2)[eset]))
  }
})

module.scores.df.cohort2 <- do.call("rbind",scores_list_cohort2)
colnames(module.scores.df.cohort2) <- c("pathway","sample","module.score","celltype")
module.scores.df.cohort2 <- merge(module.scores.df.cohort2,fgsea_res_cohort2[,c("pathway","celltype","size","leadingEdge")],by=c("pathway","celltype"),all.x=T)

# add meta data and transform the dataframe into esetlist
g_cohort2.idcelltype <- readRDS("input/SchulteSchrepping/cohort2/cluster_labels_res.0.4.rds")

severity_gsva_esetlist_cohort2 <- df_to_gsva_elist(gsva_df = module.scores.df.cohort2, meta_eset = g_cohort2.idcelltype)
saveRDS(severity_gsva_esetlist_cohort2, file.path(OUT_DIR, "severe-mild_module_score_gsva_filtered_samples_genes_cohort2.rds"))

