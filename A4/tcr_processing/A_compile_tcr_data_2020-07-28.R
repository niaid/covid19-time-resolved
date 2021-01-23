library(tidyverse)

#These can be found on GEO GSE161918
SEURAT_IN_PATH <- "input/brescia_paper1_seurat.rds"

#These can be found on GEO GSE161918
B1_TENX_IN_PATH <- "input/B1_TCR_filtered_contig_annotations.csv"
B2_TENX_IN_PATH <- "input/B2_TCR_filtered_contig_annotations.csv"
B3_TENX_IN_PATH <- "input/B3_TCR_filtered_contig_annotations.csv"

OUT_PATH <- "data/CITE5p/all_batches/tcr/2020_07_28_tenx_filtered_anno_wfiltered_celltype.rds"
dir.create(dirname(OUT_PATH))

obj <- readRDS(SEURAT_META_IN_PATH)
meta <- obj@meta.data

b1_tenx <- read_csv(B1_TENX_IN_PATH)
b2_tenx <- read_csv(B2_TENX_IN_PATH)
b3_tenx <- read_csv(B3_TENX_IN_PATH)

batches <- unique(meta$Batch)
names(batches) <- batches

tenx_list <- list(B1 = b1_tenx, B2 = b2_tenx, B3 = b3_tenx)

tenx_list <- lapply(names(tenx_list), function(batch){
  tenx <- tenx_list[[batch]]
  if(batch %in% c("B2", "B3")){
    barcode_pre <- substr(tenx$barcode, 1, 16)
    barcode_suffix <- sapply(strsplit(tenx$barcode, "-"), `[[`, 2)
    for(i in seq_along(barcode_suffix)){
      if(nchar(barcode_suffix[[i]]) == 1){
        barcode_suffix[[i]] <- paste0("0", barcode_suffix[[i]])
      }
    }
    barcodes <- paste(barcode_pre, barcode_suffix, sep = "-")
  }else{
    barcodes <- tenx$barcode
  }
  tenx$tenx_barcode <- barcodes
  tenx$barcodeBatch2 <- paste(barcodes, gsub("B", "", batch), sep = "_")
  return(tenx)
})
names(tenx_list) <- batches

prefixes <- sapply(strsplit(meta$barcodeBatch, split = "_"), `[[`, 1)
stopifnot(all(prefixes == meta$NewBarcode))

meta <- meta %>% mutate(barcodeBatch2 = paste(NewBarcode, gsub("B", "", Batch), sep = "_"))

tenx_combined <- bind_rows(tenx_list)

out_dat <- left_join(tenx_combined, meta, by = "barcodeBatch2")

out_dat <- out %>% select(-barcodeBatch2)

saveRDS(out_dat, OUT_PATH)

