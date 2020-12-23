library(limma)

FIT_IN_PATH <- snakemake@input[[1]]

CONTRAST_MAT_IN_PATH <- snakemake@input[[2]]

CONTRAST_FIT_OUT_PATH <- snakemake@output[[1]]

fit <- readRDS(FIT_IN_PATH)

contrast_mat <- readRDS(CONTRAST_MAT_IN_PATH)

cfit <- contrasts.fit(fit, contrasts = contrast_mat)

cfit <- eBayes(cfit)

saveRDS(cfit, CONTRAST_FIT_OUT_PATH)
