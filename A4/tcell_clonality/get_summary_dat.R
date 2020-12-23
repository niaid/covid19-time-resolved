library(tidyverse)
library(lmerTest)

FIT_LIST_IN_PATH <- snakemake@input[[1]]

SUMMARY_DAT_OUT_PATH <- snakemake@output[[1]]

fit_list <- readRDS(FIT_LIST_IN_PATH)

summary_dat_list <- lapply(fit_list, function(fit){
  as.data.frame(summary(fit)$coefficients) %>%
    rownames_to_column(var = "coeff") %>%
    rename(pval = `Pr(>|t|)`, t = `t value`)

})

combined_dat <- bind_rows(summary_dat_list, .id = "measure")

write_tsv(combined_dat, SUMMARY_DAT_OUT_PATH)
