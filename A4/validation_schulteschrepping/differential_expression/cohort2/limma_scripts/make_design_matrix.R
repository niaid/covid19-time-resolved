library(limma)
library(edgeR)

str(snakemake)

snakemake@source("../utility_functions/fix_single_level_contrasts.R")

#For testing
#DGELIST_IN_PATH <- "data/CITE5p/B1/time_DE/pseudobulk_dgelists_normalized/B_Naive.lambdalight-pseudobulk_dgelist_normalized.rds"
#FORMULA <- ~ days_since_hospitalized + age + sex


DGELIST_IN_PATH <- snakemake@input[[1]]

DESIGN_OUT_PATH <- snakemake@output[[1]]

FORMULA <- snakemake@params[["formula"]]
cat("Formula = ")
print(FORMULA)
FORMULA <- as.formula(FORMULA)

dge <- readRDS(DGELIST_IN_PATH)

# Make design matrix -----------------------------------------------
meta <- dge$samples
print(colnames(meta))

#drop unused factor levels
for(nm in colnames(meta)){
  if(is.factor(meta[[nm]])){
    meta[[nm]] <- droplevels(meta[[nm]])
  }
}

FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)

design <- model.matrix(FORMULA, meta)

colnames(design) <- gsub(":", "X", colnames(design))

saveRDS(design, DESIGN_OUT_PATH)
