library(Biobase)

ESET_IN_PATH <- snakemake@input[[1]]
ESET_OUT_PATH <- snakemake@output[[1]]

KEEP_TIMEPOINTS <- snakemake@params[["keep_timepoints"]]

eset <- readRDS(ESET_IN_PATH)

eset <- eset[, eset$Timepoint %in% KEEP_TIMEPOINTS]

saveRDS(eset, ESET_OUT_PATH)
