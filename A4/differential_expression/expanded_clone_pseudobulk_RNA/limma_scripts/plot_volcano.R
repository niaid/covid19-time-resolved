library(ggplot2)
library(ggrepel)
library(readr)
#library(Cairo)

TOPTAB_IN_PATH <- snakemake@input[[1]]
FIG_OUT_PATH <- snakemake@output[[1]]

FORMULA <- snakemake@params[["formula"]]
COEF <- snakemake@wildcards[["coef"]]
CELLTYPE <- snakemake@wildcards[["celltype"]]

plot_title <- paste0(CELLTYPE, "\n",
                     FORMULA, "\n",
                     "Coef: ", COEF)

toptab <- read_tsv(TOPTAB_IN_PATH)

p = ggplot(data = toptab, aes(x = logFC, y= -log(P.Value))) +
  geom_point(aes(color=adj.P.Val < .05, alpha = 0.9), show.legend = FALSE) +
  geom_text_repel(data = toptab[rank(toptab$adj.P.Val) < 40, ], aes(label = gene)) +
  scale_color_manual(values=c("darkgrey", "firebrick4")) +
  theme_classic() +
  ggtitle(plot_title)

pdf(FIG_OUT_PATH)
print(p)
dev.off()
