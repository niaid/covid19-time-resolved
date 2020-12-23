
### CITEseq cell frequency stat 
### untility funs ############################################################################
### adding basic meta data
addcellmeta <- function(df){
  ### read in metadata
  load("input/covid19.metadata.paper1.RData")
  covid19.samples$sample_id <- paste(covid19.samples$subject_id, covid19.samples$visit, sep = "_")
  covid19.samples.pbmc <- covid19.samples[covid19.samples$material_type == "PBMC",]
  covid19.samples.pbmc$severity_outcome <- paste(covid19.samples.pbmc$severity,
                                                 covid19.samples.pbmc$outcome,
                                                 sep = "_")
  PC1 <- readRDS("input/citeseq.patient.end.points.RDS")
  PC1$PC1class <-  ifelse(PC1$PC1 > median(PC1$PC1), "PC1_high", "PC1_low")
  
  df$Timepoint <- sapply(str_split(df$Var1, pattern = "_"), function(x)x[3])
  df$Timepoint <- replace(df$Timepoint, df$Timepoint == "PBMC", "HC")
  df$sample_id <- sapply(str_split(df$Var1, pattern = "_", 2), function(x)x[2])
  df$Batch <- sapply(str_split(df$Var1, pattern = "_"), function(x)x[1])
  df$Subject <- sapply(str_split(df$Var1, pattern = "_"), function(x)x[2])
  
  df$severity <- plyr::mapvalues(x = df$sample_id, 
                                 from = covid19.samples.pbmc$sample_id, to = covid19.samples.pbmc$severity)
  df$severity <- as.character(df$severity)
  df$severity <- replace(df$severity, df$Timepoint == "HC", "HC")
  df$severity <- factor(df$severity, levels = c("HC","Moderate", "Severe", "Critical"))
  df$outcome <- plyr::mapvalues(x = df$sample_id, 
                                from = covid19.samples.pbmc$sample_id, to = covid19.samples.pbmc$outcome)
  df$outcome <- as.character(df$outcome)
  df$outcome <- replace(df$outcome, df$Timepoint == "HC", "HC")
  df$severity_outcome <- paste(df$severity, df$outcome, sep = "_")
  
  df$days_since_symptoms_onset <- plyr::mapvalues(x = df$sample_id, 
                                                  from = covid19.samples.pbmc$sample_id, 
                                                  to = covid19.samples.pbmc$days_between_sample_drawn_and_symptom_onset)
  df$days_since_symptoms_onset <- as.character(df$days_since_symptoms_onset)
  df$days_since_symptoms_onset <- as.numeric(replace(df$days_since_symptoms_onset, 
                                                     df$Timepoint == "HC", 0))
  df$PC1 <- plyr::mapvalues(x = df$Subject, 
                            from = PC1$subject_id, to = PC1$PC1)
  df$PC1 <- as.character(df$PC1)
  df$PC1 <- replace(df$PC1, df$Timepoint == "HC", NA)
  df$PC1 <- as.numeric(df$PC1)
  df$PC1class <- plyr::mapvalues(x = df$Subject, from = PC1$subject_id, to = PC1$PC1class)
  df$PC1class <- replace(df$PC1class, df$Timepoint == "HC", "HC")
  df$PC1class <- factor(df$PC1class, levels = c("HC","PC1_low","PC1_high",NA))
  return(df)
}


plot5cat <- function(df, name, ratioto, width = 26, height = 18){
  # visit timepoint
  p1 <- ggplot(df, aes(x = Timepoint, y = value, color = Timepoint))+
    geom_boxplot(outlier.shape=NA)+geom_point(aes(shape = severity_outcome))+
    facet_wrap(~variable, scales = "free")+
    theme(axis.text.x = element_text(angle = 90))+
    theme_bw()
  ggsave(filename = paste(name, "samplevsTP", ratioto, "pdf", sep = "."), 
         plot = p1, device = "pdf", width = width+1.6, height = height)
  
  # severity
  p2 <- ggplot(df, aes(x = severity, y = value, color = severity))+
    geom_boxplot(outlier.shape=NA)+geom_point(aes(shape = severity_outcome))+
    facet_wrap(~variable, scales = "free")+
    theme(axis.text.x = element_text(angle = 90))+
    theme_bw()
  ggsave(filename = paste(name, "samplevsseverity", ratioto, "pdf", sep = "."), 
         plot = p2, device = "pdf", width = width+1.6, height = height)
  
  # timecourse
  tmp <- filter(df, Timepoint != "HC", PC1class %in% c("PC1_low","PC1_high"))
  tmp.HC <- filter(df, Timepoint == "HC")
  means <- plyr::ddply(tmp.HC, .(variable), summarise, median = median(value, na.rm = TRUE), 
                       quantile.25 = quantile(value, .25, na.rm = TRUE), 
                       quantile.75 = quantile(value, .75, na.rm = TRUE))
  p3 <- ggplot(tmp)+
    geom_rect(aes(xmin = 17, xmax = 23, ymin = -Inf, ymax = Inf), fill = "#EADCFA", alpha = 0.1)+
    geom_point(alpha=.4, shape=21,
               aes(x = days_since_symptoms_onset, y = value, fill=PC1class), size = 1.2)+
    scale_fill_manual(name="Severity", values = c("deepskyblue1","red")) +
    geom_line(aes(x = days_since_symptoms_onset, y = value, group = Subject), alpha = 0.08)+
    stat_summary(aes(x = days_since_symptoms_onset, y = value, group = PC1class), fun = median, alpha=0)+
    # scale_shape_manual(values = c(15:18))+
    # stat_smooth(aes(x = days_since_symptoms_onset, y = value, group = PC1class, color=PC1class), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size = 0.5)+# group = 1 overwrite group id
    stat_smooth(aes(x = days_since_symptoms_onset, y = value, group = PC1class, color=PC1class), se = FALSE, size = 0.7)+# group = 1 overwrite group id
    scale_color_manual(name="PC1 Class",values=c("deepskyblue1","red"))+
    facet_wrap(~variable, scales = "free")+
    # geom_rect(data = means, aes(xmin = -Inf, xmax = Inf, ymin = quantile.25, ymax = quantile.75), fill = "#00BA38", alpha = 0.1)+
    geom_hline(aes(yintercept = median), data = means, alpha = 0.7, linetype = "dashed", color = "#00BA38")+
    # geom_vline(xintercept = 20, color = "grey60", linetype='dashed')+
    theme(axis.text.x = element_text(angle = 90))+
    theme_bw()+
    ggtitle("green line:HC median; shaded: HC 25-75% interval")
  ggsave(filename = paste(name, "samplevsonset", ratioto, "pdf", sep = "."),
         plot = p3, device = "pdf", width = width, height = height)
  
  # PC1class
  tmp <- filter(df, !is.na(PC1class))
  tmp <- filter(tmp, Timepoint %in% c("HC", "T0"))
  p4 <- ggplot(tmp, aes(x = PC1class, y = value, color = PC1class))+
    scale_color_manual(name="PC1 Class",values=c("#00BA38", "#619CFF", "#F8766D"))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(aes(shape = severity_outcome, group = 1), position = position_jitterdodge(jitter.width = 1))+
    scale_shape_manual(values = c(15:16,3,17:18))+
    facet_wrap(~variable, scales = "free_y")+
    theme(axis.text.x = element_text(angle = 90))+
    theme_bw()
  ggsave(filename = paste(name, "samplevsPC1class", ratioto, "pdf", sep = "."),
         plot = p4, device = "pdf", width = width+1.6, height = height)
  
  # PC1 correlation
  tmp <- filter(df, Timepoint == "T0")
  p5 <- ggplot(tmp, aes(x = PC1, y = value, color = PC1))+
    geom_point(aes(shape = severity_outcome))+
    scale_shape_manual(values = c(15:18))+
    facet_wrap(~variable, scales = "free_y")+
    theme(axis.text.x = element_text(angle = 90))+
    geom_smooth(method = lm, se=FALSE, linetype = "dashed")+
    theme_bw()
  ggsave(filename = paste(name, "samplevsPC1", ratioto, "pdf", sep = "."),
         plot = p5, device = "pdf", width = width+1.6, height = height)
  
  p = list(p1,p2,p3,p4,p5)
  return(p)
}


