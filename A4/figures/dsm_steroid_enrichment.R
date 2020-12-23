library(tidyverse)
library(viridis)
library(fgsea)
library(cowplot)

FIG_fgsea_OUT_PATH <- "plots/CITE5p/all_batches/paper_figures/steroid_dsm_fgsea_curve.pdf"

gluc_model_genes <- readRDS("genesets/glucocorticoid_sets/glucResponseGenesforModel.Rds")

geneset.list <- list(glucResponseGenesforModel = gluc_model_genes)
files_toptab <- c(Mono_Classical =  "data/CITE5p/all_batches/differential_expression/2020_08_09/sample_groups/t0_covid_only/results/PC1/toptables/PC1/Unsorted-WCTcoursecelltype/Mono_Classical--model@PC1--coef@PC1--toptab.tsv", 
           NK_CD16hi = "data/CITE5p/all_batches/differential_expression/2020_08_09/sample_groups/t0_covid_only/results/PC1/toptables/PC1/Unsorted-WCTcoursecelltype/NK_CD16hi--model@PC1--coef@PC1--toptab.tsv"
           )

toptab_list <- lapply(files_toptab, function(path){
  read_tsv(path)
})

fgsea_list <- lapply(toptab_list, function(toptab){
  t.stat <- toptab$t
  names(t.stat) <- toptab$gene
  
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=15,
                    maxSize=500,
                    nperm=100000)

})


plot_list <- list()
for(nm in names(toptab_list)){
  pval <- fgsea_list[[nm]] %>% pull(pval) %>%
          round(3)
  ranks <- toptab_list[[nm]] %>%
          select(gene, t) %>%
          deframe()
  plot_list[[nm]] <- 
          plotEnrichment(stats = ranks, pathway = gluc_model_genes) +
          ylim(c(-.6, .6)) +
          ggtitle(paste(nm, "\nDSM\ngluccocorticoid signature")) +
           annotate("text",  x=Inf, y = Inf, label = paste("p =", pval), vjust=1, hjust=1)
}

#combined_dat %>% filter(pathway == path & coef == coeff.) %>% as.data.frame()
p <- plot_grid(plotlist = plot_list, ncol = 2)

ggsave(plot = p, filename = FIG_fgsea_OUT_PATH, height = 4, width = 8)
