library(edgeR)
library(tidyverse)
library(Biobase)

#IN_PATH <- "data/CITE5p/all_batches/differential_expression/2020_07_20/t0_covid_only/pseudobulk_dgelists_unfiltered/Sorted-specific_gating.rds"
#IN_PATH <- "data/CITE5p/all_batches/differential_expression/2020_07_20/t0_covid_only/pseudobulk_dgelists/Sorted-specific_gating"

#FIG_OUT_PATH <- "~/sg_data/users/rachmaninoffn2/scratch/pbulk_lib_stats.pdf"

IN_PATH <- snakemake@input[[1]]

FIG_OUT_PATH <- snakemake@output[[1]]

MIN_CELLS_PER_SAMPLE <- snakemake@params[["min_cells_per_pool"]]

#If in IN_PATH is directory, read in contents to list, otherwise read in contents
if(file_test("-d", IN_PATH)){
  files <- list.files(IN_PATH, full.names = TRUE)

  pseudobulk_list <- lapply(files, readRDS)
  names(pseudobulk_list) <- gsub("-pseudobulk_eset.rds", "", basename(files))
}else if(file_test("-f", IN_PATH)){
  pseudobulk_list <- readRDS(IN_PATH)
}


meta_list <- lapply(pseudobulk_list, pData)

meta_dat <- bind_rows(meta_list, .id = "celltype")

meta_dat <- meta_dat %>% 
        select(n_barcodes, celltype, severity.outcome, Class, Gender) %>%
        mutate(severity.outcome = as.character(severity.outcome)) %>%
        mutate(samp_group = replace(severity.outcome, Class == "HC", "HC"))

meta_dat <- meta_dat %>%
        filter(!celltype %in% "ungated")

pdf(FIG_OUT_PATH)
p1 <- ggplot(meta_dat, aes(x = celltype, fill = samp_group)) +
        geom_bar() +
        ggtitle("N samples per celltype") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p1)

p1_gender <- ggplot(meta_dat, aes(x = celltype, fill = Gender)) +
        geom_bar() +
        ggtitle("N samples per celltype") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p1_gender)

p2 <- ggplot(meta_dat, aes(x = celltype, y = log2(n_barcodes))) + 
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = samp_group, shape = Gender)) +
        ggtitle("Number of cells per sample") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        geom_hline(yintercept = log2(MIN_CELLS_PER_SAMPLE), color = "red")
print(p2)

dev.off()
