snakemake@source("../utility_functions/limma_DE_functions.R")
# to test
#DGELIST_IN_PATH <- "data/CITE5p/B1/time_DE/pseudobulk_dgelists/B_Naive-pseudobulk_dgelist.rds"
DGELIST_IN_PATH <- snakemake@input[[1]]

DGELIST_OUT_PATH <- snakemake@output[[1]]

dge <- readRDS(DGELIST_IN_PATH)

dge <- NormalizePseudobulk(dge)

saveRDS(dge, DGELIST_OUT_PATH)
