library(tidyverse)
library(readxl)
library(viridis)

IN_DIR <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/"

FIG_OUT_PATH <-"data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/clonality_model_summary_plot.pdf"

files <- list.files(IN_DIR, full.names = F, recursive = TRUE) 
files <- grep("summary", files, value = TRUE)


files_dat <- data.frame(filename = files, 
                       sample_group = sapply(strsplit(files, split = "\\/"), `[[`, 1),
                       model = sapply(strsplit(files, split = "\\/"), `[[`, 3),
                       celltype = sapply(strsplit(files, split = "\\/"), `[[`, 5),
                       stringsAsFactors = FALSE) %>%
                       mutate(celltype = gsub("_summary_dat\\.tsv", "", celltype))


names(files) <- files
dat_list <- lapply(files, function(path){
  path <- file.path(IN_DIR, path)                       
  read_tsv(path)
})

combined_dat <- bind_rows(dat_list, .id = "filename")

combined_dat <- left_join(combined_dat, files_dat)

combined_dat <- combined_dat %>%
        mutate(mod_coeff = paste(model, coeff))

keep_comparisons <- c("days_onset days_since_onset", "PC1 PC1", 
                      "PC1_cat PC1_catPC1_low", "healthy_vs_covid ClassHC")

combined_dat <- combined_dat %>%
        filter(mod_coeff %in% keep_comparisons)

#p <- ggplot(combined_dat, aes(x = measure, y = mod_coeff)) +
#        geom_point(aes(size = -log10(pval), color = t, shape = pval < .05)) +
#        scale_color_viridis() +
#        facet_wrap(celltype~sample_group, scales = "free") +
#        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p <- ggplot(combined_dat, aes(x = celltype, y = mod_coeff)) +
        geom_point(aes(size = -log10(pval), color = t, shape = pval < .05)) +
        scale_color_viridis() +
        facet_wrap(~measure, scales = "free") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p, filename = FIG_OUT_PATH, width = 16, height = 16)


