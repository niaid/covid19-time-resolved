library(tidyverse)

TENX_IN_PATH <- snakemake@input[[1]]
TENX_IN_PATH <- "data/CITE5p/all_batches/tcr/2020_07_28_tenx_filtered_anno_wfiltered_celltype.rds"

OUT_PATH <- snakemake@output[[1]]

KEEP_CELLTYPES <- snakemake@params[["keep_celltypes"]]
KEEP_BATCHES <- snakemake@params[["keep_batches"]]
KEEP_SORTED <- snakemake@params[["keep_sorted"]]

tenx <- readRDS(TENX_IN_PATH)

clones <- tenx %>% 
        filter(productive, high_confidence, is_cell) %>%
        select(barcodeBatch, Donor, Timepoint, Sample, Batch, WCTcoursecelltype,
               raw_clonotype_id, Sorted, Class, Age, Gender, severity.outcome, severity) %>% distinct()

clones <- clones %>% 
        filter(!is.na(WCTcoursecelltype)) %>%
        filter(WCTcoursecelltype %in% KEEP_CELLTYPES) %>%
        filter(Batch %in% KEEP_BATCHES) %>%
        filter(Sorted == KEEP_SORTED) 

saveRDS(clones, OUT_PATH)


