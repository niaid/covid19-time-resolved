library(tidyverse)

snakemake@source("util/diversity_metrics.R")

CLONES_IN_PATH <- snakemake@input[[1]]

DIVERSITY_OUT_PATH <- snakemake@output[[1]]

MIN_CELLS <- snakemake@params[["n_cell_cutoff"]]


clones <- readRDS(CLONES_IN_PATH)

stopifnot(sum(duplicated(clones$barcodeBatch)) == 0)
clonotypes_list <- with(clones, split(raw_clonotype_id, paste(Donor, Timepoint, Batch, sep ="-")))

# subsample and calculate diversity metrics
keep_selection <- sapply(clonotypes_list, length) > MIN_CELLS
print(paste("Removing samples:", paste(names(clonotypes_list[!keep_selection]), collapse = ", ")))
clonotypes_list <- clonotypes_list[keep_selection]

set.seed(1)
subsample_size <- MIN_CELLS

diversity_results <- lapply(clonotypes_list, function(clone_vec){
  replicate(1000, diversity_all(x = sample(clone_vec, size = subsample_size), x_type = "clonotypes"))
})

diversity_dat <- diversity_results %>%
        lapply(t) %>%
        lapply(as.data.frame) %>%
        bind_rows(.id = "Sample")

diversity_dat_summary <- 
        diversity_dat %>%
        gather(key = measure, value = value, -c(Sample)) %>%
        group_by(Sample, measure) %>%
        summarise(median1000 = median(value),
                  mean1000 = mean(value),
                  sd1000 = sd(value))

diversity_dat_summary <- diversity_dat_summary %>%
        separate(col = Sample, into = c("Donor", "Timepoint", "Batch"))

saveRDS(diversity_dat_summary, DIVERSITY_OUT_PATH)

