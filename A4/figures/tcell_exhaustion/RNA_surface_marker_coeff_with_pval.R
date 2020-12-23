library(tidyverse)
library(viridis)

FIG_TOPTAB_OUT_PATH <- "plots/CITE5p/all_batches/tcr/exhaustion/FIG4S.2020_11_21_exhaustion_markers_RNA_coefficients.pdf"

markers <- c(PD1_CD279 = "PDCD1", TIM3_CD366 = "HAVCR2", LAG3_CD223 = "LAG3", 
              TIGIT = "TIGIT", CTLA4_CD152 = "CTLA4",BTLA_CD272 = "BTLA", 
              `2B4_CD244` ="CD244", CD39 ="ENTPD1")

IN_DIR <- "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_07_31/sample_groups"

files <- list.files(IN_DIR, recursive = TRUE)
files <- grep("toptab", files, value = TRUE)
files <- c(grep("PC1", files, value = TRUE), grep("healthy_vs_covid", files, value = TRUE))

toptab_list <- lapply(files, function(path){
  read_tsv(file.path(IN_DIR, path))
})
#names(toptab_list) <- sapply(strsplit(files, "\\/"), `[[`, 7)
names(toptab_list) <- files

#lapply(toptab_list, function(dat){
#  grep("ENTPD1", dat$gene, value = T)
#})

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
        filter(coef %in% c("ClassCOVID", "PC1_catPC1_high", "PC1"))
combined_dat <- combined_dat %>%
        mutate(celltype2 = paste(gsub("-WCTcoursecelltype", "", cell_anno), celltype))

combined_dat <- combined_dat %>%
        filter(grepl("CD8", celltype2))

combined_dat <- combined_dat %>%
        mutate(feature = names(markers)[match(gene, markers)])

combined_dat <- combined_dat %>%
        filter(coef %in% c("ClassCOVID", "PC1_catPC1_high")) %>%
        filter(celltype == "CD8_Mem_expandedTRUE") %>% 
        filter(startsWith(cell_anno, "Sorted"))

combined_dat <- combined_dat %>%
        mutate(feature = names(markers)[match(gene, markers)])

p <- ggplot(combined_dat , aes(x = feature)) +
        ylab("logFC") +
        geom_col(aes(y = logFC, fill = factor(sign(logFC)))) +
        geom_text(aes(y = 1, color = P.Value < .05, label = paste("p =", round(P.Value, 3))), size = 3)+
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        facet_wrap(~coef, ncol = 1) +
        #scale_color_viridis() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position = "none")+
              #axis.line = element_line(colour = 'black', size = 2)) + 
        ylim(c(-1.1, 1.1)) +
        ggtitle("CD8 Mem: RNA Pseudobulk DE exhaustion markers ") +
        theme(plot.title = element_text(face = "bold", size = 8))

ggsave(plot = p, FIG_TOPTAB_OUT_PATH, height = 5, width = 3)
