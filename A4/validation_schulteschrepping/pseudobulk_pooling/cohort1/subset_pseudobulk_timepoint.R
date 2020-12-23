print(str(snakemake))

PBULKLIST_IN_PATH <- snakemake@input[[1]]

PBULKLIST_OUT_PATH <- snakemake@output[[1]]

KEEP_TIMEPOINTS <- snakemake@params[["keep_timepoints"]]
print("keep timepoints")
print(KEEP_TIMEPOINTS)

pbulklist <- readRDS(PBULKLIST_IN_PATH)

pbulklist <- lapply(pbulklist, function(dge){
  dge[, dge$samples$Timepoint %in% KEEP_TIMEPOINTS]
})

saveRDS(pbulklist, PBULKLIST_OUT_PATH)

