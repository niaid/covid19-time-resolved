library(tidyverse)

IN_PATH <- "data/CITE5p/all_batches/tcr/2020_07_28_tenx_filtered_anno_wfiltered_celltype.rds"

OUT_PATH <- "data/CITE5p/all_batches/tcr/2020_07_28_tcell_expansion_df.tsv"

tenx <- readRDS(IN_PATH)

out_dat <- tenx %>%
        filter(high_confidence, is_cell, productive, !is.na(Donor)) %>%
        select(barcodeBatch, Donor, Timepoint, Batch, Sorted, raw_clonotype_id) %>%
        distinct()
        #select(c(colnames(tenx)[1:27], colnames(tenx)[58: ncol(tenx)]))

out_dat <- out_dat %>%
        mutate(sample_clone = paste(Donor,Timepoint, Batch, Sorted, raw_clonotype_id, sep = "-"))

clone_counts <- table(sample_clone = out_dat$sample_clone)

clone_counts <- as.data.frame(clone_counts, stringsAsFactors = FALSE) %>%
        mutate(expanded = Freq > 1)


out_dat <- out_dat %>% left_join(clone_counts)

write_tsv(out_dat, OUT_PATH)

