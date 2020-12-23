library(tidyverse)
library(viridis)

IN_DIR <- "data/CITE5p/all_batches/differential_expression_cell_freq/2020_07_17"

FIG_OUT_PATH <- "plots/CITE5p/all_batches/differential_expression_cell_freq/cellf_freq_limma_results_interaction_compare_2020_07_19.pdf"
dir.create(dirname(FIG_OUT_PATH))

sample_groups <- list.files(IN_DIR)
names(sample_groups) <- sample_groups


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

combined_dat <- bind_rows(toptab_list, .id = "sample_group")

combined_dat <- combined_dat %>%
        mutate(ID=gsub("-model", "", ID)) %>%
        mutate(ID=gsub("-coef", "", ID)) %>%
        separate(ID, sep = "@", into = c("total_or_parent", "model", "coef"))%>%
        mutate(mod_coef = paste(model, coef, sep = "_")) %>%
        rename(cell = gene)

combined_dat <- combined_dat %>%
        filter(cell != "ungated")

combined_dat <- combined_dat %>%
        filter(grepl("PC", coef) | grepl("PLS", coef))

p <- ggplot(combined_dat, aes(x = coef, y = cell)) +
        geom_point(aes(size = -log10(adj.P.Val), color = t, shape = adj.P.Val < .05)) +
        scale_color_gradientn(colors= c("green","blue", "white", "red", "yellow"),
                              values = c(0, .4, .5, .6, 1), limits = c(-8, 8),
                              breaks = seq(-8,8, by =2))+
        #scale_colour_viridis_c() +
        facet_wrap(~total_or_parent, scales = "free") +
        #theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p, filename =FIG_OUT_PATH, height = 8, width =14)
