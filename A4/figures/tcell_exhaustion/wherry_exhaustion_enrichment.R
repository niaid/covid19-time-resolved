library(tidyverse)

genesets <- readRDS("genesets/kegg_go_btm_reactome_foointerferon.rds")

FIG_OUT_PATH <- "plots/CITE5p/all_batches/tcr/exhaustion/2020_09_01_wherry_exhaustion_enrichment.pdf"

wherry_exhaust_up <- genesets$wherry_2007_exhaustion_up
wherry_exhaust_dn <- genesets$wherry_2007_exhaustion_dn

genesets_sub <- genesets[setdiff(names(genesets), c("wherry_2007_exhaustion_dn", "wherry_2007_exhaustion_up"))]

genesets <- genesets[sapply(genesets, length) > 0]

universe <- unique(unlist(genesets))


source("/hpcdata/sg/sg_data/users/rachmaninoffn2/dev/R/module_fisher.R")

res_wherry_up <- genesetFET(hits = wherry_exhaust_up, background = universe, all.sets = genesets)
res_wherry_dn <- genesetFET(hits = wherry_exhaust_dn, background = universe, all.sets = genesets)


pdf(FIG_OUT_PATH, height = 2, width = 7)
p1 <- ggplot(res_wherry_dn %>% filter(rank(p.adj) < 10 & !grepl("wherry", geneset_name)), 
             aes(y = geneset_name, x = -log10(p.adj) ))+
        geom_col() +
        ggtitle("Wherry 2007 exhaustion\nDown")
print(p1)

p1 <- ggplot(res_wherry_up %>% filter(rank(p.adj) < 10 & !grepl("wherry", geneset_name)), 
             aes(y = geneset_name, x = -log10(p.adj) ))+
        geom_col() +
        ggtitle("Wherry 2007 exhaustion\nUp")
print(p1)

dev.off()
#file.remove(FIG_OUT_PATH)

#res_wherry_up %>% filter(p.adj < .05) %>% select(-overlap)
#res_wherry_dn %>% filter(p.adj < .05) %>% select(-overlap)
