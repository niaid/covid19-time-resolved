library(tidyverse)
library(ggrepel)

IN_DIR <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/all_timepoints_covid_only/dat_list"

FIG_OUT_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/related_plots/pc1_group_days_onset_vs_clonality_shownames.pdf"
dir.create(dirname(FIG_OUT_PATH))

files <- list.files(IN_DIR)
names(files) <- gsub("_diversity.rds", "", files)

dat_ll <- lapply(files, function(path){
  path <- file.path(IN_DIR, path)
  print(path)
  readRDS(path) %>% bind_rows(.id = "measure")
})

combined_dat <- bind_rows(dat_ll, .id = "celltype")

pdf(FIG_OUT_PATH, height =20, width = 20)

p1 <- combined_dat %>% 
        ggplot(aes(x = days_since_onset, y = median1000)) +
        geom_point() +
        #geom_text(aes(color = severity.outcome, label = Donor)) +
        geom_line(aes(group = Donor), alpha = .2) +
        geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
        facet_wrap(celltype~measure, scales = "free", nrow = 4)
print(p1)

p2 <- combined_dat %>% filter(Timepoint %in% c("HC", "T0")) %>%
        ggplot(aes(x = PC1_cat, y = median1000)) +
        geom_boxplot(outlier.shape = NA) +
        #geom_jitter(height = 0, aes(color = severity.outcome)) +
        geom_text_repel(aes(label = gsub("HGR000", "", Donor))) +
        facet_wrap(celltype~measure, scales = "free")
print(p2)

dev.off()
keep_donors <- c("HGR0000069", "HGR0000042", "HGR0000101", "HGR0000094")
combined_dat %>% filter(Donor %in% keep_donors) %>%
        select(Donor, Gender, Age, severity, outcome, PC1_cat) %>%
        distinct()



