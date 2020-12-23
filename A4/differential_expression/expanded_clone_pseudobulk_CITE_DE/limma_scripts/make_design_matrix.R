library(limma)
library(Biobase)

str(snakemake)

#For testing
#FORMULA <- ~ days_since_hospitalized + age + sex
snakemake@source("../utility_functions/fix_single_level_contrasts.R")


ESET_IN_PATH <- snakemake@input[[1]]

DESIGN_OUT_PATH <- snakemake@output[[1]]

FORMULA <- snakemake@params[["formula"]]
cat("Formula = ")
print(FORMULA)
FORMULA <- as.formula(FORMULA)

eset <- readRDS(ESET_IN_PATH)

# Make design matrix -----------------------------------------------
meta <- pData(eset)
print(colnames(meta))

#drop unused factor levels
for(nm in colnames(meta)){
  print(nm)
  if(is.factor(meta[[nm]])){
    meta[[nm]] <- droplevels(meta[[nm]])
  }
}

FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
design <- model.matrix(FORMULA, meta)

colnames(design) <- gsub(":", "X", colnames(design))

saveRDS(design, DESIGN_OUT_PATH)
