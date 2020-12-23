
fit <- readRDS("data/CITE5p/all_batches/differential_expression_cell_freq/2020_07_24/t0_covid_only/results/PC1_cat/limma_voom_fit/WCT_merged_pcnt_parent-model@PC1_cat-fit.rds")

fit$t["CD8_Mem_EM.TE", ]
fit$p.value["CD8_Mem_EM.TE", ]
