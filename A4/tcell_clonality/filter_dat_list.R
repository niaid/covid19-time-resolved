#This removes samples with NA feaures from dat_list

DAT_LIST_IN_PATH <-  "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/t0_plus_healthy/results/healthy_vs_covid/filtered_dat_list/CD4_Mem-Sorted_diversity.rds"
DAT_LIST_IN_PATH <- snakemake@input[[1]]
DAT_LIST_OUT_PATH <- snakemake@output[[1]]

FORMULA <- snakemake@params[["formula"]]

FORMULA <- as.formula(FORMULA)
FORMULA <- median1000~ Class + Age  + (1|batch)

terms_obj <- terms(FORMULA)

terms_vec <- rownames(attr(terms_obj, "factors"))
fixed_eff_terms <- terms_vec[!grepl("\\|", terms_vec)]

dat_list <- readRDS(DAT_LIST_IN_PATH)

dat_list <- lapply(dat_list, function(dat){
  keep_rows <- complete.cases(dat[, fixed_eff_terms])
  dat[keep_rows, ]
})

saveRDS(dat_list, DAT_LIST_OUT_PATH)
