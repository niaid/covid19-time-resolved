
IN_PATH <- "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/sample_groups/t0_plus_healthy/results/healthy_vs_covid/filtered_dat_list/CD4_Mem-Sorted_diversity.rds"

dat_list <- readRDS(IN_PATH)

lapply(dat_list, `[[`, "Class")
lapply(dat_list, `[[`, "Age")
lapply(dat_list, `[[`, "batch")
