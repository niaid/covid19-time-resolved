library(tidyverse)
library(ggrepel)

IN_DIR <- "data/CITE5p/all_batches/differential_expression_cell_freq/2020_08_25/"
FIG_OUT_PATH <- "plots/CITE5p/all_batches/paper_figures/FIG4/2020_08_31/FIG4_cd4_8_volcano.pdf"
FIG_OUT_PATH2 <- "plots/CITE5p/all_batches/paper_figures/FIG4/2020_08_31/FIG4_cd4_8_volcano_sep_cd4_8.pdf"
dir.create(dirname(FIG_OUT_PATH), recursive = TRUE)


files <- list.files(IN_DIR, full.names = F, recursive = TRUE) 
files <- grep("toptab", files, value = TRUE)

#files <- files[sapply(strsplit(files, "\\/"), length) == 7]

files_dat <- data.frame(filename = files, 
                       sample_group = sapply(strsplit(files, split = "\\/"), `[[`, 1),
                       model = sapply(strsplit(files, split = "\\/"), `[[`, 3),
                       coefficient = sapply(strsplit(files, split = "\\/"), `[[`, 5),
                       #cell_anno = sapply(strsplit(files, split = "\\/"), `[[`, 6),
                       celltype_mod_coef = sapply(strsplit(files, split = "\\/"), `[[`, 6),
                       stringsAsFactors = FALSE) %>%
                       mutate(celltype = sapply(strsplit(celltype_mod_coef, split = "--"), `[[`, 1))
          
keep_celltypes <- "TotalCD4and8_LE_genes_20200824"
keep_coefs <- c("COVID-Healthy", "PC1", "PC1high-low")

files_dat <- files_dat %>%
        filter(celltype %in% keep_celltypes & coefficient %in% keep_coefs)

dat_list <- lapply(files_dat$filename, function(path){
  path <- file.path(IN_DIR, path)                       
  read_tsv(path)
})
names(dat_list) <- files_dat$filename
dat_list <- lapply(dat_list, function(dat){
  dat %>% filter(!gene %in% c("TotalCD4", "TotalCD8")) %>%
          mutate(celltype2 = ifelse(grepl("CD4", gene), "CD4", "CD8"))
})

pdf(FIG_OUT_PATH, height = 4, width = 7)
for(i in seq_along(dat_list)){
        toptab <- dat_list[[i]]
        plot_title <- paste(files_dat$sample_group[[i]], files_dat$coefficient[[i]])
p = ggplot(data = toptab, aes(x = logFC, y= -log10(P.Value))) +
  geom_point(aes(color=adj.P.Val < .05, alpha = 0.9), show.legend = FALSE) +
  #geom_text_repel(data = toptab[rank(toptab$adj.P.Val) < 40, ], aes(label = gene)) +
  geom_text_repel(aes(label = gene), size = 2) +
  scale_color_manual(values=c("darkgrey", "firebrick4")) +
  geom_hline(yintercept = -log10(.05)) +
  geom_vline(xintercept = 0) +
  ylim(c(0, max(-log10(toptab$P.Value) + .5))) +
  xlim(c(-max(abs(toptab$logFC)), max(abs(toptab$logFC)))) +
  theme_classic() +
  facet_wrap(~celltype2) + 
  ggtitle(plot_title)

print(p)

}

dev.off()
#file.remove(FIG_OUT_PATH)

#plot them separately

pdf(FIG_OUT_PATH2, height = 1.7, width = 1.7)
#for(i in seq_along(dat_list)){
for(i in c(3)){
        toptab <- dat_list[[i]]

        ylims <- c(0, max(-log10(toptab$P.Value) + .5))
        xlims <- c(-max(abs(toptab$logFC)), max(abs(toptab$logFC)))
        for(cell in c("CD4", "CD8")){
          #plot_title <- paste(cell, files_dat$sample_group[[i]], files_dat$coefficient[[i]])
          plot_title <- paste(cell, files_dat$coefficient[[i]], sep = "\n")

          p = ggplot(data = toptab %>% filter(celltype2 == cell), aes(x = logFC, y= -log10(P.Value))) +
                geom_point(aes(color=adj.P.Val < .05, alpha = 0.9), show.legend = FALSE) +
                #geom_text_repel(data = toptab[rank(toptab$adj.P.Val) < 40, ], aes(label = gene)) +
                geom_text_repel(aes(label = gene), size = 1.5) +
                #geom_text_repel(aes(label = gene)) +
                scale_color_manual(values=c("darkgrey", "firebrick4")) +
                geom_hline(yintercept = -log10(.05)) +
                geom_vline(xintercept = 0) +
                ylim(ylims) +
                xlim(xlims) +
                theme_classic(base_size = 6) +
                #theme(base_size = 2)
                ggtitle(plot_title)
          print(p)
          }

}
dev.off()
#file.remove(FIG_OUT_PATH)
