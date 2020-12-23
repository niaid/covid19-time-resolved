library(tidyverse)

source("scripts/nr/CITE5p/all_batches/paper1_figures/color_schemes.R")

CD8MEM_SORTED_IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/diversity_metrics/CD8_Mem-Sorted_diversity.rds"

META_IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_meta.tsv"



meta <- read_tsv(META_IN_PATH)

diversity_dat <- readRDS(CD8MEM_SORTED_IN_PATH)

combined_dat <- left_join(diversity_dat, meta)

combined_dat <- combined_dat %>%
        mutate(PC1_cat = as.character(PC1_cat)) %>%
        mutate(PC1_cat = replace(PC1_cat, Class == "HC", "HC")) %>%
        mutate(PC1_cat= factor(PC1_cat, levels = c("HC", "PC1_low", "PC1_high"))) %>%
        filter(!is.na(PC1_cat)) %>%
        filter(measure == "simpson") %>%
        filter(!grepl("CHI", Donor))

#cd8mem_sorted_datlist <- readRDS("data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/t0_covid_only/dat_list/CD8_Mem-Sorted_diversity.rds")
#
#cd8mem_simpson <- cd8mem_sorted_datlist$simpson
#
#x1 <- 
#        combined_dat %>% 
#        mutate(dtb = paste(Donor, Timepoint, Batch)) %>%
#        select(dtb, median1000)
#
#x2 <-
#        cd8mem_simpson %>% 
#        mutate(dtb = paste(Donor, Timepoint, Batch)) %>%
#        select(dtb, median1000)
#
#x <- full_join(x1, x2, by = "dtb")


FIG_OUT_PATH <- "plots/CITE5p/all_batches/paper_figures/FIG4/2020_09_04/fig4_clonality_boxplots.pdf"
dir.create(dirname(FIG_OUT_PATH))

pdf(FIG_OUT_PATH, height =2.5, width = 4)

p1 <- combined_dat %>% 
        ggplot(aes(x = days_since_onset, y = median1000)) +
        geom_point(aes(color = severity.outcome), size = 4) +
        scale_color_manual(values = severity.color) +
        geom_smooth(method = lm, formula = y ~ splines::bs(x, 3)) +
        ylab("Clonality (Simpson Index)") +
        theme_bw() +
        ggtitle("Clonality in Sorted CD8 Mem")
print(p1)

p2 <- combined_dat %>% 
        filter(Timepoint %in% c("HC", "T0")) %>%
        ggplot(aes(x = PC1_cat, y = median1000, color = PC1_cat)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(height = 0, aes(shape = severity.outcome2), size = 3) +
        scale_shape_manual(values = severity.shape) +
        scale_color_manual(values = PC1class.color) +
        theme_bw() +
        ylab("Clonality (Simpson Index)") +
        ggtitle("Clonality in Sorted CD8 Mem")
print(p2)

dev.off()

