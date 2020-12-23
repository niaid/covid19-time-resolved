library(fgsea)
library(tidyverse)

DE_RES_IN_PATH <- "data/externalData/SchulteSchrepping/differential_expression_cohort1_no_cutoffs/sample_groups/t0_plus_healthy/results/healthy_vs_covid/toptables/COVID-Healthy/id.celltype/8_pDCs--model@healthy_vs_covid--coef@COVID-Healthy--toptab.tsv"

LE_IN_PATH <- "genesets/pdc_apoptosis_le/pDC_PC1_covid_LE.rds"

FIG_OUT_PATH <- "plots/validations/SchulteSchrepping/pdc_apoptosis_no_cell_number_cutoff/pdc_apoptosis_cohort1_using_our_LE_no_cutoff_covid-healthy.pdf"
dir.create(dirname(FIG_OUT_PATH), recursive = TRUE)

toptab <- read_tsv(DE_RES_IN_PATH)
genes <- readRDS(LE_IN_PATH)

stat_vec <- toptab %>% select(gene, t) %>% deframe()

enrich_res <- fgsea(pathways = list(pdc_apoptosis_le = genes),
                    stats = stat_vec, nperm = 10000)

pval <- round(enrich_res$pval, 5)

p <- plotEnrichment(pathway = genes, stats = stat_vec) +
             annotate("text",  x=Inf, y = Inf, label = paste("p =", pval), vjust=1, hjust=1) +
             ggtitle("cohort1 pDC apoptosis using LE from our data - no cell # cutoff")
ggsave(plot = p, filename = FIG_OUT_PATH, width = 4, height = 4)
