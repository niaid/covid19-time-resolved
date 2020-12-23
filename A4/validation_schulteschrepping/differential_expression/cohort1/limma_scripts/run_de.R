snakemake@source("../utility_functions/limma_DE_functions.R")
library(limma)
library(edgeR)

#DGELIST_IN_PATH <- "data/CITE5p/B1/time_DE/pseudobulk_dgelists_normalized/B_Naive.lambdalight-pseudobulk_dgelist_normalized.rds"
#DESIGN_IN_PATH <- "data/CITE5p/B1/time_DE/onset/pseudobulk/design_matrices/B_Naive.lambdalight-onset-design.rds"
#Run differential expression

str(snakemake)

DGELIST_IN_PATH <- snakemake@input[["dgelist"]]
DESIGN_IN_PATH <- snakemake@input[["design"]]

PATIENT_ID_COL <- snakemake@params[["patient_id_col"]]

FIT_OUT_PATH <- snakemake@output[[1]]

dgelist <- readRDS(DGELIST_IN_PATH)
stopifnot(!is.null(dgelist$samples$norm.factors))

design <- readRDS(DESIGN_IN_PATH)

#model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
dgelist <- dgelist[, rownames(design)]

fit <- RunVoomLimma(dgelist, design_matrix = design, add_random_effect = TRUE, 
                    block_this_variable = PATIENT_ID_COL,
                    do_contrast_fit = FALSE, my_contrast_matrix = NULL)

saveRDS(fit, FIT_OUT_PATH)
