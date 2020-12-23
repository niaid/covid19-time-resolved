library(tidyverse)
library(ggpubr)
library(Biobase)

eset <- readRDS("data/CITE5p/all_batches/differential_expression_cell_freq/2020_08_25/t0_plus_healthy/subsetted_expressionsets/TotalCD4and8_LE_genes_20200824-eset.rds")


mat <- exprs(eset)

dat <- as.data.frame(t(mat))

pdat <- pData(eset)


pdat <- pdat %>%
        mutate(batch_donor_time = paste(Batch, Donor, Timepoint, sep = "_")) %>%
        mutate(PC1_cat = as.character(PC1_cat)) %>%
        mutate(PC1_cat = replace(PC1_cat, Class == "HC", "HC")) %>%
        mutate(PC1_cat = factor(PC1_cat, levels = c("HC", "PC1_low", "PC1_high")))


dat <- dat %>%
        rownames_to_column(var = "batch_donor_time") %>%
        gather(key = celltype, value = value, -batch_donor_time)

dat <- dat %>% left_join(pdat)

p <- ggplot(dat %>% filter(!is.na(PC1_cat)), aes(x = PC1_cat, y = value)) +
        geom_boxplot(aes(color = PC1_cat), outlier.shape = NA) +
        geom_jitter(height = 0) +
        facet_wrap(~celltype, scales = "free_y") +
        stat_compare_means(comparisons = list(c("PC1_high", "PC1_low")))

FIG_OUT_PATH <- "plots/scratch/check_cellfreq_pc1_res.pdf"
dir.create(dirname(FIG_OUT_PATH))
ggsave(plot = p, filename = FIG_OUT_PATH, height = 20, width = 20)

mtx_orig <- readRDS("data/CITE5p/all_batches/cell_freqs/TotalCD4and8.LG.20200824.merge_cell_UnSort_WCTmerged_to_parent_mtx.Rds")

mtx <- mtx_orig %>% mutate(batch_donor_time = paste(Batch, Subject, Timepoint, sep = "_")) %>%
        `rownames<-`(.$batch_donor_time)


mtx_pc1_cat <- mtx %>% select(batch_donor_time, PC1class, PC1) %>%
        rename(PC1_mtx = PC1)

eset_pc1_cat <- pData(eset) %>% select(batch_subj_time, PC1_cat, PC1) %>% 
        rename(batch_donor_time = batch_subj_time, PC1_eset = PC1)

combined_dat <- full_join(mtx_pc1_cat, eset_pc1_cat) %>% 
        mutate(same_pc1_class = as.character(PC1class) == gsub("_", "", as.character(PC1_cat))) %>%
        mutate(same_pc1= round(PC1_mtx, 6) == round(PC1_eset, 6))
combined_dat


combined_dat %>% filter(!same_pc1_class)


mtx <- mtx[, startsWith(colnames(mtx), "Total")]

mtx <- t(as.matrix(mtx))

mtx <- mtx[rownames(mat), colnames(mat)]

identical(mat, mtx)

library(txtplot)

txtplot(mat, mtx)
mat == mtx
(as.vector(mat), as.vector(mtx))

mtx2 <- readRDS("data/CITE5p/all_batches/cell_freqs/TotalCD4and8.LG.20200824.merge_cell_UnSort_WCTmerged_to_parent_mtx (1).Rds")

mtx2$TotalCD4.LG_0
mtx_orig$TotalCD4.LG_0
identical(mtx2, mtx_orig)

mtx2$TotalCD4.LG_14
mtx_orig$TotalCD4.LG_14

