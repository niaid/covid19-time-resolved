
fit <- readRDS("data/CITE5p/all_batches/differential_expression_cell_freq/2020_08_25/t0_covid_only/results/PC1_group/limma_voom_fit/TotalCD4and8_LE_genes_20200824-model@PC1_group-fit.rds")

keep_celltypes <- c("TotalCD4.LG_14", "TotalCD8.LG_14")
fit$t[keep_celltypes, ]
fit$p.value[keep_celltypes, ]

fit$t["TotalCD8.LG_0", ]

library(tidyverse)
toptab <- read_tsv("data/CITE5p/all_batches/differential_expression_cell_freq/2020_08_25/t0_covid_only/results/PC1_group/toptables/PC1high-low/TotalCD4and8_LE_genes_20200824--model@PC1_group--coef@PC1high-low--toptab.tsv")

as.data.frame(toptab)
