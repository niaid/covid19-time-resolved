qibrary(tidyverse)
library(viridis)
library(fgsea)
library(cowplot)

FIG_fgsea_OUT_PATH <- "plots/CITE5p/all_batches/tcr/exhaustion/2020_08_17_exhaustion_signature_enrichments.pdf"

GENESETS_IN_PATH <- "genesets/kegg_go_btm_reactome_foointerferon.rds"

genesets <- readRDS(GENESETS_IN_PATH)

IN_DIR <- "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups/"

files <- list.files(IN_DIR, recursive = TRUE)
files <- grep("fgsea_table", files, value = TRUE)
files <- c(grep("PC1", files, value = TRUE), grep("healthy_vs_covid", files, value = TRUE))

fgsea_list <- lapply(files, function(path){
  path <- file.path(IN_DIR, path)
  read_tsv(path,
           col_types = cols(
                       pathway = col_character(),
                       pval = col_double(),
                       padj = col_double(),
                       ES = col_double(),
                       NES = col_double(),
                       nMoreExtreme = col_double(),
                       size = col_double(),
                       leadingEdge = col_character()
                       ))
})
#names(fgsea_list) <- sapply(strsplit(files, "\\/"), `[[`, 7)
names(fgsea_list) <- files

combined_dat <- bind_rows(fgsea_list, .id = "filepath")

combined_dat <- combined_dat %>% filter(grepl("wherry", pathway))
combined_dat <- combined_dat %>% 
        mutate(cell_anno = sapply(strsplit(filepath, "\\/"), `[[`, 6)) %>%
        mutate(filename = sapply(strsplit(filepath, "\\/"), `[[`, 7)) %>%
        mutate(filename = gsub("-coef", "", filename)) %>%
        mutate(filename = gsub("-model", "", filename)) %>%
        mutate(filename = gsub("-fgsea\\.tsv", "", filename)) %>%
        separate(filename, sep = "@", into = c("celltype", "model", "coef"))

combined_dat <- combined_dat %>%
        filter(coef %in% c("ClassCOVID", "PC1_catPC1_high", "PC1")) %>%
        filter(celltype == "CD8_Mem_expandedTRUE") %>% 
        filter(startsWith(cell_anno, "Sorted"))


IN_TOPTAB_DIR <- "data/CITE5p/all_batches/expanded_tcell_pbulk_CITE_DE/2020_07_31/sample_groups"

files_toptab <- list.files(IN_TOPTAB_DIR, recursive = TRUE)
files_toptab <- grep("toptab", files_toptab, value = TRUE)
files_toptab <- c(grep("PC1", files_toptab, value = TRUE), grep("healthy_vs_covid", files_toptab, value = TRUE))

toptab_list <- lapply(files_toptab, function(path){
  read_tsv(file.path(IN_DIR, path))
})
#names(toptab_list) <- sapply(strsplit(files_toptab, "\\/"), `[[`, 7)
names(toptab_list) <- files_toptab

combined_dat_toptab <- bind_rows(toptab_list, .id = "filepath")


combined_dat_toptab <- combined_dat_toptab %>% 
        mutate(cell_anno = sapply(strsplit(filepath, "\\/"), `[[`, 6)) %>%
        mutate(filename = sapply(strsplit(filepath, "\\/"), `[[`, 7)) %>%
        mutate(filename = gsub("-coef", "", filename)) %>%
        mutate(filename = gsub("-model", "", filename)) %>%
        mutate(filename = gsub("-toptab\\.tsv", "", filename)) %>%
        separate(filename, sep = "@", into = c("celltype", "model", "coef"))

combined_dat_toptab <- combined_dat_toptab %>%
        filter(coef %in% c("ClassCOVID", "PC1", "PC1_catPC1_high")) %>%
        filter(celltype == "CD8_Mem_expandedTRUE") %>% 
        filter(startsWith(cell_anno, "Sorted"))


plot_list <- list()
for(coeff. in c("ClassCOVID", "PC1_catPC1_high")){
  for(path in unique(combined_dat$pathway)){
    pval <- combined_dat %>% filter(pathway == path & coef == coeff.) %>% pull(pval) %>%
            round(3)
    ranks <- combined_dat_toptab %>%
            filter(coef == coeff.) %>%
            select(gene, t) %>%
            deframe()
    plot_list[[paste(coeff., path, sep = "_")]] <- 
            plotEnrichment(stats = ranks, pathway = genesets[[path]]) +
            ylim(c(-.6, .6)) +
            ggtitle(paste(coeff., path, sep = "\n")) +
             annotate("text",  x=Inf, y = Inf, label = paste("p =", pval), vjust=1, hjust=1)
  }

}
#combined_dat %>% filter(pathway == path & coef == coeff.) %>% as.data.frame()
p <- plot_grid(plotlist = plot_list, ncol = 2)

ggsave(plot = p, filename = FIG_fgsea_OUT_PATH, width = 8)
