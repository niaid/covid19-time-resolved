# covid_paper1_code_nr

All code used to create figures is found in /figures

All other directories contain code that preprocesses/analyses data for figure generation.

There are few different analyses, but the three main workflows for running the preprocessing code are as follows:

tcr_processing --> pseudobulk_pooling --> differential_expression --> figures
tcr_processing --> tcell_clonality--> figures
validation_schulteschrepping/pseudobulk_pooling --> validation_schulteschrepping/differential_expression --> figures

Run scripts from each directory in the order shown above.

The analysis pipelines in tcr_processing, pseudobulk_pooling, and differential_expression, use snakemake extensively. To run this you must have both python and snakemake installed. It has only been tested with snakemake/5.8.2-Python-3.6.7; however other snakemake and python versions may still work.

To run the snakemake workflows, there are several bash scripts with the following pattern: sm_call*
These are configured to run the workflows by submitting jobs with the univa grid engine. It is possible to use another job scheduler, but this would require changing the sm_call* script and cluster_config.json to accomodate the other scheduler.  
A description of cluster submission with snakemake can be found here https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html 

Alternatively, the workflows can be run on a single machine, using the command 'snakemake -s name__of__snakefile.py' after changing to the correct directory.
For the differential_expression workflows, it is suggested to use the --keepgoing flag as some celltypes don't have sufficient samples to be run and the jobs will always fail.


See overall_workflow for an overview of the order in which to run each script in order to produce the figures.


A brief description of which script corresponds to each figure is below:

Figure 3: Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity.				
A	COVIDvsHC bubble	figures/gsea_bubble_plots.R	
B	DSM bubble	figures/gsea_bubble_plots.R	
G	corr hm and scatters	figures/fig3_heatmap_and_scatterplots.R	
				
Figure 5: Single Cell and Clonal Expansion Analysis in CD4+ and CD8+ T cells and Exhaustion Assessment in Clonal CD8+ T cells.				
C	CD8 cluster FC	figures/fig4_plot_cd4_8_fig4_cluster_volcano.R	
E	clonality DSM	figures/clonality/plot_pcgroup_onset_assoc_for_paper.R	
F	CD8 exhaustion marker FC	figures/tcell_exhaustion/surface_marker_coefficients_with_pval.R	
G	exhaustion gene GSEA curve	figures/tcell_exhaustion/exhaustion_fgsea_plots.R	
				
				
Figure S3:  Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity				
A	DSM-low TSO bubble	figures/gsea_bubble_plots.R	
B	DSM-high TSO bubble	figures/gsea_bubble_plots.R	
G	German data pDC apopt GSEA curve	figures/validations/pdc_apoptosis/*	
				
Figure S4:  Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity				
G	German data NK FA GSEA curve	figures/validations/NK_fattyacid_gsea_curves_cohorts_1_2.R	
				
Figure S5:  Exogeneous corticosteroid effect				
A	glucocorticoid GSEA curve	dsm_steroid_enrichment.R	
				
Figure S6: CD4 and CD8 exhaustion validation				
C	CD4 cluster FC	figures/fig4_plot_cd4_8_fig4_cluster_volcano.R	
F	exhaustion cutoff boxplot	figures/tcell_exhaustion/cutoff_markers_grid.R	
G	exhaustion LE enrichment	figures/tcell_exhaustion/wherry_exhaustion_enrichment.R	
H	German data exhaustion bubble	figures/validations/exhaustion_plots.R	
I	German data exhaustion GSEA curve	figures/validations/exhaustion_plots.R	
				
Figure S6: CD4 and CD8 exhaustion validation				
A	figures/gsea_bubble_plots.R	
				
Figure supplemental note: celltype specific gene expression profile				
G	figures/gsea_bubble_plots.R	
				
				
				
				

