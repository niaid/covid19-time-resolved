library(lmerTest)
print(str(snakemake))

DAT_LIST_IN_PATH <- snakemake@input[[1]]

FIT_LIST_OUT_PATH <- snakemake@output[[1]]

FORMULA <- snakemake@params[["formula"]]
Y_COL <- snakemake@params[["y_col"]]

FORMULA <- paste0(Y_COL, FORMULA)
print(FORMULA)

FORMULA <- as.formula(FORMULA)

dat_list <- readRDS(DAT_LIST_IN_PATH)

fit_list <- lapply(dat_list, function(dat){
  lmer(formula = FORMULA, data = dat)
})

saveRDS(fit_list, FIT_LIST_OUT_PATH)
