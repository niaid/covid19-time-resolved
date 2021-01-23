# A4 code

# Environment

R-3.6.0 

python-3.6.7

see environment.txt for full list of R packages needed to be installed

# Download data from GEO

Necessary data are found in GEO, GSE161918.

Download brescia_paper1_seurat.rds,  and B*_TCR_filtered_contig_annotations.csv files and put these in the input/ directory

Other required data can be found in input folder

All intermediate files necessary to produce figures will be saved in data/ folder after processing scripts are run

# How to run

See overall_workflow.sh for an overview of the order in which to run each script in order to produce the figures.

## Description

All code used to create figures is found in /figures

All other directories contain code that preprocesses/analyses data for figure generation.

There are few different analyses, but the three main workflows for running the preprocessing code are as follows:

tcr_processing --> pseudobulk_pooling --> differential_expression --> figures

tcr_processing --> tcell_clonality--> figures

validation_schulteschrepping/pseudobulk_pooling --> validation_schulteschrepping/differential_expression --> figures

Run scripts from each directory in the order shown above.

## Using Snakemake

The analysis pipelines in tcr_processing, pseudobulk_pooling, and differential_expression, use snakemake extensively. To run this you must have both python and snakemake installed. It has only been tested with snakemake/5.8.2-Python-3.6.7; however other snakemake and python versions may still work.

Rather than directly invoking snakemake with the snakemake command, there are several bash scripts with the following pattern: sm_call*

For example to run the pooling workflow for the main cell clusters, you can call sm_call_pooling with the following command

```
bash pseudobulk_pooling/pool_expanded_clone_CITE/sm_call_pooling
```

These are configured to run the workflows by submitting jobs with the univa grid engine. It is possible to use another job scheduler, but this would require changing the sm_call* script and cluster_config.json to accomodate the other scheduler.  
A description of cluster submission with snakemake can be found here https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html 

Alternatively, the workflows can be run on a single machine, by deleting the --cluster flag and arguments from sm_call* or using the command 'snakemake -s name__of__snakefile.py' after changing to the correct directory.
For the differential_expression workflows, it is suggested to use the --keepgoing flag as some celltypes don't have sufficient samples for differential expression to be run and the jobs will fail.


# Guide to figures
A brief description of which script corresponds to each figure is below:

## Figure 3: Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity.				

### A	COVIDvsHC bubble	

figures/gsea_bubble_plots.Rmd	

### B	DSM bubble	

figures/gsea_bubble_plots.Rmd	

### G	corr hm and scatters	

figures/fig3_heatmap_and_scatterplots.Rmd	
				

## Figure 5: Single Cell and Clonal Expansion Analysis in CD4+ and CD8+ T cells and Exhaustion Assessment in Clonal CD8+ T cells.				

### C	CD8 cluster FC	

figures/fig4_plot_cd4_8_fig4_cluster_volcano.Rmd	

### E	clonality DSM	

figures/clonality/plot_pcgroup_onset_assoc_for_paper.Rmd	

### F	CD8 exhaustion marker FC	

figures/tcell_exhaustion/surface_marker_coefficients_with_pval.Rmd	

### G	exhaustion gene GSEA curve	

figures/tcell_exhaustion/exhaustion_fgsea_plots.Rmd	
				
				
## Figure S3:  Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity				

### A	DSM-low TSO bubble	

figures/gsea_bubble_plots.Rmd	

### B	DSM-high TSO bubble	

figures/gsea_bubble_plots.Rmd	

### G	German data pDC apopt GSEA curve	

figures/validations/pdc_apoptosis/*	
				

## Figure S4:  Cell-type-specific Gene Expression Signatures of COVID-19 Disease Status and Severity				

### G	German data NK FA GSEA curve	

figures/validations/NK_fattyacid_gsea_curves_cohorts_1_2.Rmd	
				

## Figure S5:  Exogeneous corticosteroid effect				

### A	glucocorticoid GSEA curve	

dsm_steroid_enrichment.Rmd	
				

## Figure S6: CD4 and CD8 exhaustion validation				

### C	CD4 cluster FC	

figures/fig4_plot_cd4_8_fig4_cluster_volcano.Rmd	

### F	exhaustion cutoff boxplot	

figures/tcell_exhaustion/cutoff_markers_grid.Rmd	

### G	exhaustion LE enrichment	

figures/tcell_exhaustion/wherry_exhaustion_enrichment.Rmd	

### H	German data exhaustion bubble	

figures/validations/exhaustion_plots.Rmd	

### I	German data exhaustion GSEA curve	

figures/validations/exhaustion_plots.Rmd	
				
## Figure S6: CD4 and CD8 exhaustion validation				

### A	

figures/gsea_bubble_plots.Rmd	
				

## Supplemental table S4 and S5: celltype specific gene expression profile				

### Excel spreadsheets	

figures/gsea_bubble_plots.Rmd	
				
				
				
				

