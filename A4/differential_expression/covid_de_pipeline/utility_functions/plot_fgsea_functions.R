#Modified from github/MattPM/de_workflow/downstream_analysis_functions.r  on July 12, 2020

# For GSEA bubble plot instead of heatmap (plots NES and P value and only + enrichment shown)
# merge GSEA list into dataframe 
RbindGseaResultList = function(gsea_result_list, NES_filter = 0, padj_filter = 0.05){
  
  score = lapply(gsea_result_list, function(x){ 
    x = x %>% 
      select(pathway, padj, NES, celltype) %>% 
      filter( NES > NES_filter ) 
  })
  score = do.call(rbind, score) %>% filter(padj < padj_filter) %>% mutate(n_logp = -log10(padj))
  return(score)
}




GSEABubblePlot2 = function(rbind_gsea_result_dataframe, save_path = NULL, main = "", width = 8.5, height = 7.2, show_plot = FALSE) { 
  p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, size = -log10(padj), fill = NES)) + 
    geom_point(shape = 21) + 
    scale_fill_viridis_c() + 
    theme_bw() +
    scale_x_discrete(position = "top") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "right") + 
    labs(size = '-log10(ajdusted P value)', fill = 'Normalized Enrichment Score') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))+
    ggtitle(main)

  if(!is.null(save_path)){
    ggsave(p, filename = paste0(save_path), width = width, height = height)
  }

  if(show_plot){
    print(p)
  }
  return(p)
  
}
