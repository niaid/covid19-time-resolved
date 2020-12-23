### utility funs ##########################################################################################
### plot time since onset geneset enrichment score
library(ggplot2)
plot_onset_module_score_line <- function(input_gsva_esetlist, LE_score_used = "GSVA PC1 LE score",
                                              fig_out_pdf = "modulescore_time_allcelltype.pdf", 
                                              celltype = celltypes,
                                              width = 8.5, height = 6.5) {
  p <- list()
  pdf(paste(FIG_OUT_PATH,fig_out_pdf,sep = ""), width = width, height = height)
  for(cell in celltype){
    if(grepl("dblt", cell, ignore.case = T) | cell == "gated_out" | cell == "Unknown"){
      next()
    }
    # print(cell)
    eset <- input_gsva_esetlist[[cell]]
    dat <- t(exprs(eset))
    dat <- as.data.frame(dat) %>% select(intersect(combined_genesets, colnames(dat)))
    if(nrow(dat) < 7){
      next()
    }
    
    dat$days_since_onset <- eset$days_since_onset
    dat$pc1_group <- eset$PC1_cat
    dat$Subject <- eset$Subject
    dat$severity_outcome <- eset$severity.outcome2
    
    dat <- dat %>% 
      gather(key = module, value = score, 
             -c(days_since_onset, pc1_group, severity_outcome, Subject)) %>%
      filter(!is.na(pc1_group))
    
    tmp.HC <- filter(dat, pc1_group == "HC" & Subject != "CHI014")
    means <- tmp.HC %>% dplyr::group_by(module) %>% 
      dplyr::summarise(median = median(score),
                       quantile.25 = quantile(score, .25, na.rm = TRUE),
                       quantile.75 = quantile(score, .75, na.rm = TRUE))
    
    tmp.covid <- filter(dat, pc1_group %in% c("PC1_low", "PC1_high"))
    severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
    p1 <- ggplot(tmp.covid, aes(x = days_since_onset, y = score)) + 
      # geom_rect(aes(xmin = 17, xmax = 23, ymin = -Inf, ymax = Inf), fill = "#EADCFA", alpha = 0.1)+
      geom_point(alpha=.4,shape=21,aes(fill=pc1_group),size=2) + 
      scale_fill_manual(name="Severity",values = c("deepskyblue1","red")) +
      geom_line(aes(group = Subject), alpha = 0.1)+
      ggtitle(paste(cell, LE_score_used, sep = " ")) + 
      scale_color_manual(name="PC1 Class",values=c("deepskyblue1","red")) +
      geom_smooth(se = T,aes(color=pc1_group),alpha=0.1) + 
      geom_hline(aes(yintercept = median), data = means, alpha = 0.7, linetype = "dashed", color = "#00BA38")+
      facet_wrap(~module, scales = "free", ncol = 4) +
      geom_vline(xintercept = 20, color = "grey60", linetype='dashed') +
      theme_bw()
    print(p1)
    p[[cell]] <- p1
  }
  dev.off()
  return(p)
}


## functions for highlight the TSO d17-23 period
plot_onset_module_score_highlight <- function(input_gsva_esetlist, LE_score_used = "GSVA PC1 LE score",
                                    fig_out_pdf = "modulescore_time_allcelltype.pdf", 
                                    celltype = celltypes,
                                    width = 8.5, height = 6.5) {
  p <- list()
  pdf(paste(FIG_OUT_PATH,fig_out_pdf,sep = ""), width = width, height = height)
  for(cell in celltype){
    if(grepl("dblt", cell, ignore.case = T) | cell == "gated_out" | cell == "Unknown"){
      next()
    }
    # print(cell)
    eset <- input_gsva_esetlist[[cell]]
    dat <- t(exprs(eset))
    dat <- as.data.frame(dat) %>% select(intersect(combined_genesets, colnames(dat)))
    if(nrow(dat) < 7){
      next()
    }
    
    dat$days_since_onset <- eset$days_since_onset
    dat$pc1_group <- eset$PC1_cat
    dat$Subject <- eset$Subject
    dat$severity_outcome <- eset$severity.outcome2
    
    dat <- dat %>% 
      gather(key = module, value = score, 
             -c(days_since_onset, pc1_group, severity_outcome, Subject)) %>%
      filter(!is.na(pc1_group))
    
    tmp.HC <- filter(dat, pc1_group == "HC" & Subject != "CHI014")
    means <- tmp.HC %>% dplyr::group_by(module) %>% 
      dplyr::summarise(median = median(score),
                       quantile.25 = quantile(score, .25, na.rm = TRUE),
                       quantile.75 = quantile(score, .75, na.rm = TRUE))
    
    tmp.covid <- filter(dat, pc1_group %in% c("PC1_low", "PC1_high"))
    severity.color <- c("Critical-Alive"="#374E55FF","Critical-Deceased"="#DF8F44FF","Moderate-Alive"="#00A1D5FF" ,"Severe-Alive"="#B24745FF","HC" = "#79AF97FF")
    p1 <- ggplot(tmp.covid, aes(x = days_since_onset, y = score)) + 
      geom_rect(aes(xmin = 17, xmax = 23, ymin = -Inf, ymax = Inf), fill = "#EADCFA", alpha = 0.1)+
      geom_point(alpha=.4,shape=21,aes(fill=pc1_group),size=2) + 
      scale_fill_manual(name="Severity",values = c("deepskyblue1","red")) +
      geom_line(aes(group = Subject), alpha = 0.1)+
      ggtitle(paste(cell, LE_score_used, sep = " ")) + 
      scale_color_manual(name="PC1 Class",values=c("deepskyblue1","red")) +
      geom_smooth(se = T,aes(color=pc1_group),alpha=0.1) + 
      geom_hline(aes(yintercept = median), data = means, alpha = 0.7, linetype = "dashed", color = "#00BA38")+
      facet_wrap(~module, scales = "free", ncol = 4) +
      # geom_vline(xintercept = 20, color = "grey60", linetype='dashed') +
      theme_bw()
    print(p1)
    p[[cell]] <- p1
  }
  dev.off()
  return(p)
}



