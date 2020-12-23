library(limma)
library(tidyverse)

FIT_IN_PATH <- snakemake@input[[1]]
TOPTAB_OUT_PATH <- snakemake@output[[1]]


CONTRAST_COL <- snakemake@wildcards[["coef"]]

fit <- readRDS(FIT_IN_PATH)

toptab <- topTable(fit, number = nrow(fit), coef = CONTRAST_COL) %>%
        rownames_to_column(var = "gene")

write_tsv(toptab, TOPTAB_OUT_PATH)

