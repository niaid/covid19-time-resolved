#### Script descriptions

## 1.preprocessData.R

	This script takes the processed data (HDF5 files) and demuxlet files from the GEO submission and creates a Seurat object, performs filtering of the cells as well as hashtag demultiplexing and DSB normalization of the CITEseq data.
	
## 2.clustering.R

	This script takes the Seurat object output by 1.preprocessData.R, and performs 2 rounds of clustering, first within the batches, then within each celltype, to obtain the final cell subsets.
	
## 3.UMAPandProteinVis.Figure2.R

	This script takes the processed Seurat object from GEO as input, and generates UMAP plots and protein heatmaps that are used in Figure 2 of the manuscript.
	
## 4.TcellReclustering.R

	This script takes the processed Seurat object from GEO as input, and generates the T cell RNA clustering analysis and visualizations that are used in Figure 5 and S6 of the manuscript.
	
## 5.externaldatasets.R

	This script takes the processed Seurat object from GEO and external datasets processed data (from Schulte-Schrepping et al. Cell. 2020 Sep 17; 182(6): 1419–1440.e23 PMID: 32810438) as input, and projects matching cluster names onto the external datasets.
	
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.2 (Maipo)

Matrix products: default
BLAS/LAPACK: 

locale:
[1] en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scico_1.1.0        DescTools_0.99.29  genefilter_1.66.0  magrittr_1.5       ggridges_0.5.1    
 [6] viridis_0.5.1      viridisLite_0.3.0  pheatmap_1.0.12    parallelDist_0.2.4 forcats_0.4.0     
[11] stringr_1.4.0      purrr_0.3.2        readr_1.3.1        tidyr_0.8.3        tibble_2.1.1      
[16] ggplot2_3.3.2      tidyverse_1.2.1    matrixStats_0.54.0 dplyr_0.8.0.1      Seurat_3.2.2      

loaded via a namespace (and not attached):
  [1] readxl_1.3.1          backports_1.1.4       plyr_1.8.4            igraph_1.2.4.1       
  [5] lazyeval_0.2.2        splines_3.6.1         listenv_0.7.0         digest_0.6.18        
  [9] htmltools_0.3.6       gdata_2.18.0          memoise_1.1.0         tensor_1.5           
 [13] cluster_2.0.8         ROCR_1.0-7            limma_3.40.6          globals_0.12.4       
 [17] annotate_1.62.0       modelr_0.1.4          RcppParallel_4.4.4    colorspace_1.4-1     
 [21] blob_1.1.1            rvest_0.3.3           ggrepel_0.8.1         haven_2.1.0          
 [25] crayon_1.3.4          RCurl_1.95-4.12       jsonlite_1.6          spatstat_1.64-1      
 [29] spatstat.data_1.4-3   survival_2.44-1.1     zoo_1.8-6             glue_1.4.2           
 [33] polyclip_1.10-0       gtable_0.3.0          leiden_0.3.1          future.apply_1.3.0   
 [37] BiocGenerics_0.30.0   abind_1.4-5           scales_1.0.0          mvtnorm_1.0-10       
 [41] DBI_1.0.0             miniUI_0.1.1.1        Rcpp_1.0.1            xtable_1.8-4         
 [45] reticulate_1.12       bit_1.1-14            rsvd_1.0.2            stats4_3.6.1         
 [49] htmlwidgets_1.3       httr_1.4.0            gplots_3.0.1.1        RColorBrewer_1.1-2   
 [53] ica_1.0-2             pkgconfig_2.0.2       XML_3.98-1.19         uwot_0.1.8           
 [57] deldir_0.1-29         tidyselect_1.1.0      labeling_0.3          rlang_0.4.7          
 [61] reshape2_1.4.3        later_0.8.0           AnnotationDbi_1.46.0  munsell_0.5.0        
 [65] cellranger_1.1.0      tools_3.6.1           cli_1.1.0             generics_0.0.2       
 [69] RSQLite_2.1.1         broom_0.5.2           yaml_2.2.0            goftest_1.2-2        
 [73] npsurv_0.4-0          bit64_0.9-7           fitdistrplus_1.0-14   caTools_1.17.1.2     
 [77] RANN_2.6.1            pbapply_1.4-2         future_1.14.0         nlme_3.1-139         
 [81] mime_0.6              xml2_1.2.0            compiler_3.6.1        rstudioapi_0.10      
 [85] plotly_4.9.0          png_0.1-7             lsei_1.2-0            spatstat.utils_1.17-0
 [89] stringi_1.4.3         RSpectra_0.15-0       lattice_0.20-38       Matrix_1.2-17        
 [93] pillar_1.3.1          lmtest_0.9-36         RcppAnnoy_0.0.13      data.table_1.12.2    
 [97] cowplot_1.0.0         bitops_1.0-6          irlba_2.3.3           httpuv_1.5.1         
[101] patchwork_0.0.1.9000  R6_2.4.0              promises_1.0.1        KernSmooth_2.23-15   
[105] gridExtra_2.3         IRanges_2.18.2        codetools_0.2-16      boot_1.3-20          
[109] MASS_7.3-51.3         gtools_3.8.1          assertthat_0.2.1      withr_2.1.2          
[113] sctransform_0.3       S4Vectors_0.22.0      expm_0.999-4          mgcv_1.8-28          
[117] parallel_3.6.1        hms_0.4.2             grid_3.6.1            rpart_4.1-15         
[121] Rtsne_0.15            Biobase_2.44.0        shiny_1.3.2           lubridate_1.7.4      
