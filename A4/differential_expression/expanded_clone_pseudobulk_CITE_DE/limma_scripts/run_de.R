snakemake@source("../utility_functions/limma_DE_functions.R")
library(limma)
library(edgeR)
library(Biobase)

#eset_IN_PATH <- "data/CITE5p/B1/time_DE/pseudobulk_esets_normalized/B_Naive.lambdalight-pseudobulk_eset_normalized.rds"
#DESIGN_IN_PATH <- "data/CITE5p/B1/time_DE/onset/pseudobulk/design_matrices/B_Naive.lambdalight-onset-design.rds"
#Run differential expression

str(snakemake)

ESET_IN_PATH <- snakemake@input[["eset"]]
DESIGN_IN_PATH <- snakemake@input[["design"]]

PATIENT_ID_COL <- snakemake@params[["patient_id_col"]]

FIT_OUT_PATH <- snakemake@output[[1]]

eset <- readRDS(ESET_IN_PATH)

design <- readRDS(DESIGN_IN_PATH)

#model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the eset
eset <- eset[, rownames(design)]

fit <- RunLimma(eset, design_matrix = design, add_random_effect = TRUE, 
                    block_this_variable = PATIENT_ID_COL,
                    do_contrast_fit = FALSE, my_contrast_matrix = NULL)

saveRDS(fit, FIT_OUT_PATH)
