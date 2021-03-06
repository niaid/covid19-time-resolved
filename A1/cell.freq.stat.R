### get the cell frequency of unsorted CITEseq data
### the input is the Seurat object
library(Seurat) #load Seurat 3.1
library(matrixStats)
library(plyr)
library(tidyverse)
library(openxlsx)
library(reshape2)
source("util_fun/cellfreq.funs.R")
### unsorted #####################################################################################
# get cell frequency and ratio of total from Seurat object
# Download the Seurat object first and put it into input folder
merge <- readRDS("input/brescia_paper1_seurat.rds")
merge.filtered <- subset(merge, subset = celltypeQC == TRUE)
merge.UnSort <- subset(merge.filtered, subset = Sorted == "N")

### read in metadata
load("input/covid19.metadata.paper1.RData")
covid19.samples$sample_id <- paste(covid19.samples$subject_id, covid19.samples$visit, sep = "_")
covid19.samples.pbmc <- covid19.samples[covid19.samples$material_type == "PBMC",]
covid19.samples.pbmc$severity_outcome <- paste(covid19.samples.pbmc$severity,
                                               covid19.samples.pbmc$outcome,
                                               sep = "_")
PC1 <- readRDS("input/citeseq.patient.end.points.RDS")
PC1$PC1class <-  ifelse(PC1$PC1 > median(PC1$PC1), "PC1_high", "PC1_low")


### cell ratios UnSort
merge_cell_UnSort_WCTmerged <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$WCTmergedcelltype)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq))
merge_cell_UnSort_WCTcourse <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$WCTcoursecelltype)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq))
merge_cell_UnSort_adjusted <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$adjustedcelltype)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq))
merge_cell_UnSort_adjcourse <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$coursecelltype)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq))

# saveRDS(merge_cell_UnSort_WCTmerged, "output/merge.UnSort.WCTmerged.sample.rds")
# saveRDS(merge_cell_UnSort_WCTcourse, "output/merge.UnSort.WCTcourse.sample.rds")
# saveRDS(merge_cell_UnSort_adjusted, "output/merge.UnSort.adjusted.sample.rds")
# saveRDS(merge_cell_UnSort_adjcourse, "output/merge.UnSort.adjustedcourse.sample.rds")

### Unsort summary -- WCTmergedcelltype (high resolution cell populations) ##################################
# parse cell frequency data
### get adjcourse file first for calculating parent ratio
merge_cell_UnSort_adjcourse <- addcellmeta(merge_cell_UnSort_adjcourse)

### WCTmergedcelltype
merge_cell_UnSort_WCTmerged <- addcellmeta(merge_cell_UnSort_WCTmerged)

### get ratio to their parent celltypes, also merge course celltype columns
### (B subsets to CD19+, T subsets to CD4/CD8+, DC/baso/others to myeloid)
merge_cell_UnSort_WCTmerged_mtx <- merge_cell_UnSort_WCTmerged %>% 
  dcast(Var1+Timepoint+sample_id+Batch+severity+
          days_since_symptoms_onset+PC1class+PC1+Subject+severity_outcome ~ Var2, value.var = "ratio")

# get freq matrix to filter out CD19, CD4, CD8, NK, Mono with freq<10 cells
merge_cell_UnSort_adjcourse_Freq_mtx <- merge_cell_UnSort_adjcourse %>% 
  dcast(Var1 ~ Var2, value.var = "Freq") %>%
  mutate(CD19 = B+PB)


merge_cell_UnSort_adjcourse_mtx <- merge_cell_UnSort_adjcourse %>% 
  dcast(Var1 ~ Var2, value.var = "ratio") %>%
  mutate(CD19 = B+PB) %>%
  left_join(merge_cell_UnSort_adjcourse_Freq_mtx[,c("Var1","CD4","CD8","CD19","NK","Mono")], by = "Var1") %>%
  mutate("CD19" = replace(.$CD19.x, .$CD19.y < 10, values = NA)) %>%
  mutate("CD4" = replace(.$CD4.x, .$CD4.y < 10, values = NA)) %>%
  mutate("CD8" = replace(.$CD8.x, .$CD8.y < 10, values = NA)) %>%
  mutate("NK" = replace(.$NK.x, .$NK.y < 10, values = NA)) %>%
  mutate("Mono" = replace(.$Mono.x, .$Mono.y < 10, values = NA))


merge_cell_UnSort_WCTmerged_mtx <- left_join(merge_cell_UnSort_WCTmerged_mtx, 
                                             merge_cell_UnSort_adjcourse_mtx[,c("Var1","CD4","CD8","CD19","NK","Mono")], 
                                             "Var1")

merge_cell_UnSort_WCTmerged_to_parent_mtx <- merge_cell_UnSort_WCTmerged_mtx %>%
  mutate_at(vars(contains('B_')), ~(.)/CD19) %>%
  mutate_at(vars(contains('CD4_')), ~(.)/CD4) %>%
  mutate_at(vars(contains('CD8_')), ~(.)/CD8) %>%
  mutate_at(vars(contains('NK_')), ~(.)/NK) %>%
  mutate_at(vars(contains('Mono_')), ~(.)/Mono) %>%
  mutate(Treg = Treg/CD4)

# get the NAs reverted back for the populations themselves Freq to total
identical(merge_cell_UnSort_adjcourse_mtx$Var1, merge_cell_UnSort_WCTmerged_mtx$Var1)
merge_cell_UnSort_WCTmerged_mtx <- merge_cell_UnSort_WCTmerged_mtx %>% 
  mutate("CD4" = merge_cell_UnSort_adjcourse_mtx$CD4.x,
         "CD8" = merge_cell_UnSort_adjcourse_mtx$CD8.x,
         "CD19" = merge_cell_UnSort_adjcourse_mtx$CD19.x,
         "NK" = merge_cell_UnSort_adjcourse_mtx$NK.x,
         "Mono" = merge_cell_UnSort_adjcourse_mtx$Mono.x)

merge_cell_UnSort_WCTmerged_to_parent <- merge_cell_UnSort_WCTmerged_to_parent_mtx %>% 
  melt(id = c("Var1", "Timepoint", "sample_id", "Batch", "severity", "days_since_symptoms_onset","PC1","PC1class","Subject","severity_outcome"))

# plot5cat(merge_cell_UnSort_WCTmerged_to_parent, "merge.UnSort.WCTmerged", "toparent")

### Unsort summary -- WCTcoursecelltype (coarser level cell populations) ######################
merge_cell_UnSort_WCTcourse <- addcellmeta(merge_cell_UnSort_WCTcourse)

### get ratio to their parent celltypes, also merge coarser celltype columns
### (B subsets to CD19+, T subsets to CD4/CD8+, DC/baso/others to myeloid)
merge_cell_UnSort_WCTcourse_mtx <- merge_cell_UnSort_WCTcourse %>% 
  dcast(Var1+Timepoint+sample_id+Batch+severity+
          days_since_symptoms_onset+PC1+PC1class+Subject+severity_outcome ~ Var2, value.var = "ratio")

merge_cell_UnSort_WCTcourse_mtx <- left_join(merge_cell_UnSort_WCTcourse_mtx, merge_cell_UnSort_adjcourse_mtx[,c("Var1","CD4","CD8","CD19","NK","Mono")], 
                                             "Var1")
merge_cell_UnSort_WCTcourse_to_parent_mtx <- merge_cell_UnSort_WCTcourse_mtx %>%
  mutate_at(vars(contains('B_')), ~(.)/CD19) %>%
  mutate_at(vars(contains('CD4_')), ~(.)/CD4) %>%
  mutate_at(vars(contains('CD8_')), ~(.)/CD8) %>%
  mutate_at(vars(contains('NK_')), ~(.)/NK) %>%
  mutate_at(vars(contains('Mono_')), ~(.)/Mono) %>%
  mutate(Treg = Treg/CD4)


# get the NAs reverted back for the populations themselves Freq to total
identical(merge_cell_UnSort_adjcourse_mtx$Var1, merge_cell_UnSort_WCTcourse_mtx$Var1)
merge_cell_UnSort_WCTcourse_mtx <- merge_cell_UnSort_WCTcourse_mtx %>% 
  mutate("CD4" = merge_cell_UnSort_adjcourse_mtx$CD4.x,
         "CD8" = merge_cell_UnSort_adjcourse_mtx$CD8.x,
         "CD19" = merge_cell_UnSort_adjcourse_mtx$CD19.x,
         "NK" = merge_cell_UnSort_adjcourse_mtx$NK.x,
         "Mono" = merge_cell_UnSort_adjcourse_mtx$Mono.x)

merge_cell_UnSort_WCTcourse_to_total <- merge_cell_UnSort_WCTcourse_mtx %>% 
  melt(id = c("Var1", "Timepoint", "sample_id", "Batch", "severity", "days_since_symptoms_onset","PC1","PC1class","Subject","severity_outcome"))
# plot5cat(merge_cell_UnSort_WCTcourse_to_total, "merge.UnSort.WCTcourse", "tototal")

### have a final full list of cell frequencies of all cell gated/clustered cell populations ###################
# get the cell population freq of specific gating
merge_cell_UnSort_specific <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$specific_gating)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq)) %>%
  dplyr::filter(Var2 %in% c("B_CD71", "Tfh", "Treg", "Plasmablast"))

merge_cell_UnSort_specific <- addcellmeta(merge_cell_UnSort_specific)

### get ratio to their parent celltypes
### (B subsets to CD19+, T subsets to CD4/CD8+)
merge_cell_UnSort_specific_mtx <- merge_cell_UnSort_specific %>% 
  dcast(Var1+Timepoint+sample_id+Batch+severity+days_since_symptoms_onset+PC1+PC1class+Subject+severity_outcome ~ Var2, value.var = "ratio")


# get freq matrix to filter out CD19, CD4, CD8, NK, Mono with freq<10 cells
merge_cell_UnSort_course <- data.frame(table(merge.UnSort$sample_id, merge.UnSort$gate_group_course)) %>% 
  group_by(Var1) %>% mutate(ratio = Freq/sum(Freq))

merge_cell_UnSort_course_Freq_mtx <- merge_cell_UnSort_course %>% 
  dcast(Var1 ~ Var2, value.var = "Freq") %>%
  mutate(CD19 = B+Plasmablast)

merge_cell_UnSort_course_mtx <- merge_cell_UnSort_course %>% 
  dcast(Var1 ~ Var2, value.var = "ratio") %>%
  mutate(CD19 = B+Plasmablast) %>%
  left_join(merge_cell_UnSort_course_Freq_mtx[,c("Var1","B","CD4","CD8","CD19","NK","Mono")], by = "Var1") %>%
  mutate("CD19" = replace(.$CD19.x, .$CD19.y < 10, values = NA)) %>%
  mutate("CD4" = replace(.$CD4.x, .$CD4.y < 10, values = NA)) %>%
  mutate("CD8" = replace(.$CD8.x, .$CD8.y < 10, values = NA)) %>%
  mutate("NK" = replace(.$NK.x, .$NK.y < 10, values = NA)) %>%
  mutate("Mono" = replace(.$Mono.x, .$Mono.y < 10, values = NA))

merge_cell_UnSort_mtx <- left_join(merge_cell_UnSort_specific_mtx, merge_cell_UnSort_course_mtx[,c("Var1","CD4","CD8","CD19","NK","Mono")], by = "Var1")

merge_cell_UnSort_to_parent_mtx <- merge_cell_UnSort_mtx %>%
  mutate_at(vars(contains('B_')), ~(.)/CD19) %>%
  mutate(Plasmablast = Plasmablast/CD19,
         Tfh = Tfh/CD4,
         Treg = Treg/CD4)

merge_cell_UnSort_to_parent_mtx_specific <- select(merge_cell_UnSort_to_parent_mtx, 
                                                   Var1,Plasmablast,B_CD71,Treg,Tfh)
colnames(merge_cell_UnSort_to_parent_mtx_specific) <- paste(colnames(merge_cell_UnSort_to_parent_mtx_specific), "gated", sep = "_")


# get final freq full list of all celltypes
WCTmerged_parent_selected <- setdiff(colnames(merge_cell_UnSort_WCTmerged_to_parent_mtx), colnames(merge_cell_UnSort_WCTcourse_to_parent_mtx))
merge_cell_UnSort_WCTmerged_to_parent_mtx_selected <- select(merge_cell_UnSort_WCTmerged_to_parent_mtx, Var1, all_of(WCTmerged_parent_selected))

WCTcousre_parent_selected <- c("B_Naive", "B_Mem", "PB_Plasmablasts", "Treg", 
                               "CD4_Mem", "CD8_Mem", "CD8_Naive", "CD4_Naive")
merge_cell_UnSort_WCTcourse_to_parent_mtx_selected <- select(merge_cell_UnSort_WCTcourse_to_parent_mtx,Var1, all_of(WCTcousre_parent_selected))

merge_cell_UnSort_to_parent_mtx_selected <- 
  left_join(merge_cell_UnSort_WCTcourse_to_parent_mtx_selected, merge_cell_UnSort_to_parent_mtx_specific, by = c("Var1" = "Var1_gated")) %>% 
  left_join(., merge_cell_UnSort_WCTmerged_to_parent_mtx_selected, by = c("Var1"))
colnames(merge_cell_UnSort_to_parent_mtx_selected) <- paste(colnames(merge_cell_UnSort_to_parent_mtx_selected), "to_parent", sep = "_")

# remove the ones from to_total already included in to_parent
WCTcourse_total_selected <- setdiff(colnames(merge_cell_UnSort_WCTcourse_mtx),
                                    WCTcousre_parent_selected)

# remove "Unknown","Dblt","dim"
WCTcourse_total_selected <- WCTcourse_total_selected[-which(WCTcourse_total_selected %in% c("Unknown","Dblt","dim"))]
merge_cell_UnSort_to_total_mtx_selected <- select(merge_cell_UnSort_WCTcourse_mtx,
                                                  all_of(WCTcourse_total_selected))

merge_cell_UnSort_mtx_selected <- left_join(merge_cell_UnSort_to_total_mtx_selected, 
                                            merge_cell_UnSort_to_parent_mtx_selected, 
                                            by = c("Var1" = "Var1_to_parent"))
write.csv(merge_cell_UnSort_mtx_selected, "output/final_full_cell_freq_UnSort_mtx.20200818.csv")





