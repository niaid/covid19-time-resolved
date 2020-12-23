library(tidyverse)
library(fgsea)
library(cowplot)

IN_DIR <- "data/externalData/SchulteSchrepping/"

GENESETS_IN_PATH <- "genesets/kegg_go_btm_reactome_foointerferon.rds"

genesets <- readRDS(GENESETS_IN_PATH)

files <- list.files(IN_DIR, recursive = TRUE)

files <- c(grep("fgsea_table", files, value = T), grep("toptables", files, value = T))

files_split <- strsplit(files, "/")

files_dat <- lapply(1:9, function(i){sapply(files_split, `[[`, i)}) %>% 
        as.data.frame(stringsAsFactors = FALSE)
colnames(files_dat) <- c("cohort", "sg", "samp_group", "rs", "model", "result_type", "coeff",
                         "annotation", "cell_mod_coeff")

files_dat <- files_dat %>% 
        mutate(full_path = file.path(IN_DIR, files)) %>%
        mutate(cell = sapply(strsplit(cell_mod_coeff, "--"), `[[`, 1)) %>%
        filter(cohort %in% c("differential_expression", "differential_expression_cohort2")) %>%
        mutate(cohort = ifelse(grepl("cohort2", cohort), "cohort2", "cohort1"))

keep_cells <- c("CD8+ T cells",
                "CD8_Mem", "12_CD8+ T cells_1",
                "13_CD8+ T cells_2", "14_CD8+ T cells_3")

files_dat <- files_dat %>%
        filter(coeff %in% c("severe-mild", "COVID-Healthy")) %>%
        filter(cell %in% keep_cells)

files_dat %>% select(annotation, cell) %>% distinct() %>% arrange(annotation)

fgsea_files_dat <- files_dat %>%
        filter(result_type == "fgsea_tables")

fgsea_dat_list <- lapply(fgsea_files_dat$full_path, read_tsv)
names(fgsea_dat_list) <- fgsea_files_dat$full_path

fgsea_dat <- bind_rows(fgsea_dat_list, .id = "full_path")

fgsea_dat <- fgsea_dat %>%
        filter(grepl("wherry", pathway)) %>%
        left_join(fgsea_files_dat)


BUBBLE_OUT_PATH <- "plots/validations/SchulteSchrepping/exhaustion/exhaustion_bubble_schulteSchrepping.pdf"
dir.create(dirname(BUBBLE_OUT_PATH))

p_bubble <- ggplot(fgsea_dat, aes(y = pathway, x = paste(cohort, annotation, cell, sep = "\n"))) +
        geom_point(aes(size = -log10(pval), color = NES, shape = pval < .05)) +
        scale_color_viridis_c() +
        facet_wrap(~coeff, scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab("") +
        theme(plot.margin= unit(c(t = 3, r = 0, b = 3, l = 0), unit = "in"))

ggsave(plot = p_bubble, filename = BUBBLE_OUT_PATH, width =10, height = 9)

LE_CSV_OUT_PATH <- "plots/validations/SchulteSchrepping/exhaustion/cohort2_leading_edge.csv"
fgsea_dat %>% filter(cohort == "cohort2") %>%
        select(pathway, coeff, annotation, cell, NES, pval, leadingEdge) %>%
        write_csv(LE_CSV_OUT_PATH)

toptab_files_dat <- files_dat %>%
        filter(result_type == "toptables")

plot_list <- list()

keep_pathways <- grep("wherry", names(genesets), value = TRUE)
keep_comparisons <- c("severe-mild", "COVID-Healthy")

OUT_DIR <- "plots/validations/SchulteSchrepping/exhaustion/2020_11_20_exhaustion_gsea_curves"
dir.create(OUT_DIR)

for(cohort_select in c("cohort1", "cohort2")){
  for(path in keep_pathways){
    for(coefficient in keep_comparisons){
      cell_dat <- fgsea_dat %>%
              filter(cohort == cohort_select) %>% 
              filter(pathway == path) %>% 
              filter(coeff == coefficient) %>% 
              select(cell, annotation) %>% unique()

        #print(1)
      for(i in seq_len(nrow(cell_dat))){
        celltype <- cell_dat[[i, "cell"]]
        anno <- cell_dat[[i, "annotation"]]

        #print(2)
        pval <- fgsea_dat %>% 
                filter(cohort == cohort_select) %>% 
                filter(pathway == path) %>% 
                filter(cell == celltype) %>% 
                filter(coeff == coefficient) %>% 
                filter(annotation == anno) %>% 
                pull(pval) %>%
                round(5)

        #print(3)
        toptab_path <- toptab_files_dat %>% 
                filter(cohort == cohort_select) %>% 
                #filter(pathway == path) %>% 
                filter(cell == celltype) %>% 
                filter(coeff == coefficient) %>% 
                filter(annotation == anno) %>% 
                pull(full_path)

        toptab <- read_tsv(toptab_path)

        ranks <- toptab %>%
                select(gene, t) %>%
                deframe()

        plot_list[[cohort_select]][[coefficient]][[path]][[paste(anno, celltype)]] <- 
                plotEnrichment(stats = ranks, pathway = genesets[[path]]) +
                ylim(c(-.6, .6)) +
                ggtitle(paste(cohort_select, coefficient, anno, celltype, path, sep = "\n")) +
                 annotate("text",  x=Inf, y = Inf, label = paste("p =", pval), vjust=1, hjust=1)
      }
      plotlist_sub <- plot_list[[cohort_select]][[coefficient]][[path]]
      fig_out_path <- paste(cohort_select, coefficient, anno, celltype, path, sep = "_")
      fig_out_path <- paste0(file.path(OUT_DIR, fig_out_path), ".pdf")
      p <- plot_grid(plotlist = plotlist_sub, nrow = 1)
      ggsave(plot = p, filename= fig_out_path, height = 4, width = 4 * length(plotlist_sub))

    }
  }
}


        fgsea_dat %>% 
                filter(cohort == cohort_select) %>% 
                filter(pathway == path) %>% 
                filter(cell == celltype) %>% 
                filter(coeff == coefficient) %>% 
                filter(annotation == anno) %>% 
                as.data.frame()

#combined_dat %>% filter(pathway == path & coef == coeff.) %>% as.data.frame()
#p <- plot_grid(plotlist = plot_list, ncol = 2)


