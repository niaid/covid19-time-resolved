print(str(snakemake))

library(limma)

DESIGN_MATRIX_PATH <- snakemake@input[[1]]

CONTRAST <- snakemake@params[["contrast"]]
CONTRAST_NAME <- snakemake@wildcards[["coef"]]

print("contrast")
print(CONTRAST)

CONTRAST_MAT_OUT_PATH <- snakemake@output[[1]]

design <- readRDS(DESIGN_MATRIX_PATH)

contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design)) 
colnames(contrast_mat) <- CONTRAST_NAME

saveRDS(contrast_mat, CONTRAST_MAT_OUT_PATH)

#tmp <- readRDS("data/CITE5p/all_batches/differential_expression/2020_08_05/sample_groups/all_timepoints_covid_only/results/onset_group/contrast_matrices/Unsorted-WCTcoursecelltype/B_Naive--model@onset_group--coef@late-early--contrast_mat.rds")
