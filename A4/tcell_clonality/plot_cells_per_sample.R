library(ggplot2)

IN_PATH <- snakemake@input[[1]]

FIG_OUT_PATH <- snakemake@output[[1]]

N_CELL_CUTOFF <- snakemake@params[["n_cell_cutoff"]]

clones <- readRDS(IN_PATH)

p <- ggplot(clones, aes(x = paste(Donor, Timepoint, Batch))) +
            geom_bar() +
            geom_hline(yintercept = N_CELL_CUTOFF, color = "red") +
            facet_wrap(~Class, scales = "free_x") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p, filename = FIG_OUT_PATH, height = 6, width = 12)
