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
### generate GSVA score files for NK.modulescore.corr.Rmd, output files also provided in input/SchulteSchrepping/cohort1
### cohort1 #######################################################################################
sample.groups <- "all_timepoints_covid_only"
DGELISTS_IN_PATH <- "input/SchulteSchrepping/cohort1/pseudobulk_dgelists_normalized/"
# get filtered data(cell #/sample and genes filtered)
files <- list.files(DGELISTS_IN_PATH, full.names = TRUE)
OUT_DIR <- "input/SchulteSchrepping/cohort1"

genesets <- readRDS("input/kegg_go_btm_reactome_foointerferon.rds")

pbulk_list <- lapply(files, readRDS)
names(pbulk_list) <- gsub("-pseudobulk_dgelist_normalized\\.rds", "", basename(files))

# list of expression data by cell type
eset_list <- lapply(pbulk_list, function(dge){
  mat <- DESeq2::varianceStabilizingTransformation(dge$counts)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})
eset_list <- unlist(eset_list)

# get fgsea results of LE genes from the severity-mild model
FGSEA_IN_PATH <- file.path("input/SchulteSchrepping/cohort1/fgsea_tables/severe-mild")
fgsea_files <- list.files(FGSEA_IN_PATH, full.names = TRUE)

fgsea_list <- lapply(fgsea_files, function(path){
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

names(fgsea_list) <- gsub("--model.+\\.tsv", "", basename(fgsea_files))
fgsea_res <- bind_rows(fgsea_list, .id = "celltype")

scores_list <- lapply(1:length(eset_list), function(eset){
  cat(names(eset_list)[eset],"\n")
  celltype.gsea.res <- subset(fgsea_res, celltype == names(eset_list)[eset])
  if (nrow(celltype.gsea.res) > 0) {
    celltype.leading.edge.genesets.df <- aggregate(leadingEdge ~ pathway,celltype.gsea.res,paste,collapse=" ")
    celltype.leading.edge.genesets <- sapply(celltype.leading.edge.genesets.df$leadingEdge,function(x){unique(unlist(strsplit(x," ")))})
    names(celltype.leading.edge.genesets) <- celltype.leading.edge.genesets.df$pathway
    module.scores <- gsva(expr = eset_list[[eset]], gset.idx.list = celltype.leading.edge.genesets, method = "gsva", parallel.sz = 16)
    return(cbind(reshape2::melt(exprs(module.scores)),celltype= names(eset_list)[eset]))
  }
})

module.scores.df <- do.call("rbind",scores_list)
colnames(module.scores.df) <- c("pathway","sample","module.score","celltype")
module.scores.df <- merge(module.scores.df,fgsea_res[,c("pathway","celltype","size","leadingEdge")],by=c("pathway","celltype"),all.x=T)

# add meta data and transform the dataframe into esetlist
g_cohort1.idcelltype <- readRDS("input/SchulteSchrepping/cohort1/id.celltype.rds")

severity_gsva_esetlist <- df_to_gsva_elist(gsva_df = module.scores.df, meta_eset = g_cohort1.idcelltype)
saveRDS(severity_gsva_esetlist, file.path(OUT_DIR, "severe-mild_module_score_gsva_filtered_samples_genes_cohort1.rds"))



