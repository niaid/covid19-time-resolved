library(tidyverse)
library(viridis)

FIG_TOPTAB_OUT_PATH <- "plots/CITE5p/all_batches/tcr/exhaustion/FIG4L.2020_11_21_exhaustion_surface_markers_CITE_pc1_healthy_coeff_values.pdf"

markers <- c(PD1_CD279 = "CD279", TIM3_CD366 = "CD366", LAG3_CD223 = "CD223", 
              TIGIT = "TIGIT", CTLA4_CD152 = "CD152",BTLA_CD272 = "CD272", 
              `2B4_CD244` ="CD244", CD39 = "CD39")

IN_DIR <- "data/CITE5p/all_batches/expanded_tcell_pbulk_CITE_DE/2020_07_31/sample_groups"

files <- list.files(IN_DIR, recursive = TRUE)
files <- grep("toptab", files, value = TRUE)
files <- c(grep("PC1", files, value = TRUE), grep("healthy_vs_covid", files, value = TRUE))

toptab_list <- lapply(files, function(path){
  read_tsv(file.path(IN_DIR, path))
})
#names(toptab_list) <- sapply(strsplit(files, "\\/"), `[[`, 7)
names(toptab_list) <- files

combined_dat <- bind_rows(toptab_list, .id = "filepath")

combined_dat <- combined_dat %>% filter(gene %in% markers)

combined_dat <- combined_dat %>% 
        mutate(cell_anno = sapply(strsplit(filepath, "\\/"), `[[`, 6)) %>%
        mutate(filename = sapply(strsplit(filepath, "\\/"), `[[`, 7)) %>%
        mutate(filename = gsub("-coef", "", filename)) %>%
        mutate(filename = gsub("-model", "", filename)) %>%
        mutate(filename = gsub("-toptab\\.tsv", "", filename)) %>%
        separate(filename, sep = "@", into = c("celltype", "model", "coef"))

combined_dat <- combined_dat %>%
        filter(coef %in% c("ClassCOVID", "PC1_catPC1_high")) %>%
        filter(celltype == "CD8_Mem_expandedTRUE") %>% 
        filter(startsWith(cell_anno, "Sorted"))

combined_dat <- combined_dat %>%
        mutate(feature = names(markers)[match(gene, markers)])

p <- ggplot(combined_dat, aes(x = feature)) +
        ylab("logFC") +
        geom_col(aes(y = logFC, fill = factor(sign(logFC)))) +
        #geom_text(aes(y = .45, color = P.Value < .05, label = paste("p =", round(P.Value, 2))), size = 1.7)+ 
        geom_text(aes(y = 0, color = P.Value < .05, label = paste("p =", round(P.Value, 2))), size = 1.7, angle = 90)+ 
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        facet_wrap(~coef, ncol = 1) +
        #scale_color_viridis() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 70, hjust=1),
              legend.position = "none")+
        theme(plot.title = element_text(size = 5, face = "bold")) +
              #axis.line = element_line(colour = 'black', size = 2)) + 
        ylim(c(-.5, .5)) +
        ggtitle("CD8 Mem: CITE Protein Pseudobulk DE exhaustion markers ")

ggsave(plot = p, FIG_TOPTAB_OUT_PATH, height = 3, width = 3)
