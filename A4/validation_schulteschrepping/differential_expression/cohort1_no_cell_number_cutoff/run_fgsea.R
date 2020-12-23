#Thsi script assumes that you have added a column called gene to the toptable instead of just using the rownames

library(limma)
library(variancePartition)
library(fgsea)
library(data.table)
#library(cowplot)

#Set paths ---------------------------------------------------------------
#Inputs
#FIT.IN.PATH <- "~/sg_data/PROJECTS/Nicaragua_Project/analysis/data/analysis_out/time_DE/sex_age_interaction/age_group_specific/dream/fit_objects/age_first_pax_1-and-2_fit.rds"
TOPTAB.IN.PATH <- snakemake@input[[1]]

COMBINED.GENESETS.IN.PATH <- snakemake@input[[2]]

#Out path
TSV.OUT.PATH <- snakemake@output[[1]]
FIG.OUT.PATH <- snakemake@output[[2]]

#Load data ---------------------------------------------------------------
geneset.list <- readRDS(COMBINED.GENESETS.IN.PATH)

toptab <- read.table(TOPTAB.IN.PATH, header = TRUE)

#Run enrichment -----------------------------------------------------------
t.stat <- toptab$t
names(t.stat) <- toptab$gene

set.seed(1)
fgseaRes <- fgsea(pathways = geneset.list, 
                  stats = t.stat,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)

topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
pdf(FIG.OUT.PATH, width = 14, height = 10)
plotGseaTable(geneset.list[topPathways], t.stat, fgseaRes, 
              gseaParam = 0.5)
dev.off()
fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
