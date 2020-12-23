library(fgsea)
library(tidyverse)

DE_RES_IN_PATH <- "data/externalData/SchulteSchrepping/differential_expression_cohort2/sample_groups/t0_plus_healthy/results/healthy_vs_covid/toptables/COVID-Healthy/cluster_labels_res.0.4/pDC--model@healthy_vs_covid--coef@COVID-Healthy--toptab.tsv"

LE_IN_PATH <- "genesets/pdc_apoptosis_le/pDC_PC1_covid_LE.rds"

FIG_OUT_PATH <- "plots/validations/SchulteSchrepping/pdc_apoptosis/pdc_apoptosis_cohort2_using_our_LE.pdf"
dir.create(dirname(FIG_OUT_PATH), recursive = TRUE)

toptab <- read_tsv(DE_RES_IN_PATH)
genes <- readRDS(LE_IN_PATH)

stat_vec <- toptab %>% select(gene, t) %>% deframe()

enrich_res <- fgsea(pathways = list(pdc_apoptosis_le = genes),
                    stats = stat_vec, nperm = 10000)

pval <- round(enrich_res$pval, 5)

p <- plotEnrichment(pathway = genes, stats = stat_vec) +
             annotate("text",  x=Inf, y = Inf, label = paste("p =", pval), vjust=1, hjust=1) +
             ggtitle("cohort2 pDC apoptosis using LE from our data")
ggsave(plot = p, filename = FIG_OUT_PATH)
