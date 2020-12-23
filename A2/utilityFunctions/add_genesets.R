library(tidyverse)
library(readxl)

geneset_list_in_path <- "/hpcdata/sg/sg_data/users/rachmaninoffn2/genesets/processed/unlisted_gene_sets.RDS"

genesets <- readRDS(geneset_list_in_path)

foo_genesets <- read_tsv("./multi_study2.txt")

hallmark_sets <- readLines("./h.all.v7.1.symbols.gmt")
hallmark_sets <- strsplit(hallmark_sets, "\t")
names(hallmark_sets) <- sapply(hallmark_sets, `[[`, 1)
hallmark_sets <- lapply(hallmark_sets, function(x){
  x[-c(1,2)]
})

foo_genesets <- as.list(foo_genesets)

ginsburg_2013_flu <- readRDS("./Ginsburg_2013_flu_tableS4.Rds")

mckinney_cd8_exhaustion <- read_tsv("./CD8_exhaustion_module.txt", col_names = FALSE)

mckinney_cd2_stim_dn <- read_excel("./41586_2015_BFnature14468_MOESM17_ESM.xlsx", sheet = 1)
mckinney_cd2_stim_up <- read_excel("./41586_2015_BFnature14468_MOESM17_ESM.xlsx", sheet = 2)

mckinney_cd2_stim_dn <- as.data.frame(mckinney_cd2_stim_dn)[-1, 1]
mckinney_cd2_stim_up <- as.data.frame(mckinney_cd2_stim_up)[-1, 1]

wherry_exhaustion_up <- read_excel("mckinney_nature_2016_14468_table_s6.xlsx", sheet =1, skip =1)
wherry_exhaustion_dn <- read_excel("./mckinney_nature_2016_14468_table_s6.xlsx", sheet =2, skip =1)

cd8_exhaustion_list <- list(mckinney_2015_cd8_aav_exhaustion_enriched = mckinney_cd8_exhaustion$X1,
                            mckinney_2015_cd8_cd2up_exhaustion_dn = mckinney_cd2_stim_up,
                            mckinney_2015_cd8_cd2dn_exhaustion_up = mckinney_cd2_stim_dn,
                            wherry_2007_exhaustion_up = toupper(wherry_exhaustion_up[[1]]),
                            wherry_2007_exhaustion_dn = toupper(wherry_exhaustion_dn[[1]]))



foo_genesets <- lapply(foo_genesets, function(geneset){
  geneset[!is.na(geneset)]
})

names(foo_genesets) <- paste0("foointerferon", names(foo_genesets))

genesets <- c(genesets,  hallmark_sets, foo_genesets, cd8_exhaustion_list, ginsburg_2013_flu)

#remove HLA genes from a few of the interferon genesets

rm_hla_genesets <- c("GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY", 
                     "GO_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA", 
                     "GO_RESPONSE_TO_INTERFERON_GAMMA")

genesets[rm_hla_genesets]

for(gset in rm_hla_genesets){
  genesets[[gset]] <- genesets[[gset]][!startsWith(genesets[[gset]], "HLA")]
}

genesets[rm_hla_genesets]
saveRDS(genesets, "kegg_go_btm_reactome_foointerferon.rds")
