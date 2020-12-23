library(tidyverse)
library(lemon)
library(viridis)

IN_DIR <- "data/CITE5p/all_batches/differential_expression_cell_freq/2020_07_24"

FIG_OUT_PATH <- "plots/CITE5p/all_batches/differential_expression_cell_freq/cell_freq_limma_results_2020_07_24.pdf"
dir.create(dirname(FIG_OUT_PATH))

sample_groups <- list.files(IN_DIR)
names(sample_groups) <- sample_groups
sample_groups <- setdiff(sample_groups, "expressionsets")


toptab_list <- lapply(sample_groups, function(sg){
  path <- file.path(IN_DIR, sg)
  files <- list.files(path, recursive = TRUE)

  files <- grep("toptab.tsv", files, value = T) 
  names(files) <- gsub("-toptab\\.tsv", "", basename(files))
  lapply(files, function(f){
    full_path <- file.path(path, f)
    read_tsv(full_path)
  }) %>% bind_rows(.id = "ID")
})
names(toptab_list) <- sample_groups

combined_dat <- bind_rows(toptab_list, .id = "sample_group")

combined_dat <- combined_dat %>%
        mutate(ID=gsub("-model", "", ID)) %>%
        mutate(ID=gsub("-coef", "", ID)) %>%
        separate(ID, sep = "@", into = c("cell_annotation", "model", "coef"))%>%
        mutate(mod_coef = paste(model, coef, sep = "_")) %>%
        rename(cell = gene)

combined_dat <- combined_dat %>%
        filter(cell != "ungated")

combined_dat <- combined_dat %>%
        mutate(pval_class = "not signif") %>%
        mutate(pval_class = replace(pval_class, P.Value < .05, "P < .05")) %>%
        mutate(pval_class = replace(pval_class, adj.P.Val < .05, "fdr < .05")) %>%
        mutate(pval_class = relevel(factor(pval_class), "not signif"))

#combined_dat %>% select(P.Value, adj.P.Val, pval_class) %>% as.data.frame() %>% head(20)

cell_annos <- unique(combined_dat$cell_annotation)
pdf(FIG_OUT_PATH, height = 16, width = 20)
for(anno in cell_annos){
        dat_sub <- combined_dat %>% filter(cell_annotation == anno)
        p <- ggplot(dat_sub, aes(x = coef, y = cell)) +
          geom_point(aes(size = -log10(adj.P.Val), color = t, shape = pval_class)) +
          scale_shape_manual(values = c(16,17, 25)) +
          scale_color_gradientn(colors= c("green","blue", "white", "red", "yellow"),
                              values = c(0, .4, .5, .6, 1), limits = c(-8, 8),
                              breaks = seq(-8,8, by =2))+
          #scale_colour_viridis_c() +
          facet_rep_wrap(~ sample_group, scales = "free_x",  repeat.tick.labels=TRUE) +
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          ggtitle(anno)
        print(p)
}

dev.off()

