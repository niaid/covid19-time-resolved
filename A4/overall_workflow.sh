#!/usr/bin/bash

# This bash function allows you to knit Rmd
function knit() {
    R -e "rmarkdown::render('$1')"
}


### Download data from GEO ###

#Download brescia_paper1_seurat.rds,  and B*_TCR_filtered_contig_annotations.csv files and put these in the input/ directory
#Other required data can be found in input folder

#Run scripts in this order

###  DE resulsts and figures for main clusters  ###

#pseudobulk_pooling 
bash pseudobulk_pooling/pool_main_clusters/sm_call_pooling

# differential expression
bash differential_expression/covid_de_pipeline/sm_call_t0_covid
bash differential_expression/covid_de_pipeline/sm_call_t0_whealthy
bash differential_expression/covid_de_pipeline/sm_call_covid_timecourse

#figures
knit figures/gsea_bubble_plots.Rmd
knit figures/fig3_heatmap_and_scatterplots.Rmd
knit figures/dsm_steroid_enrichment.Rmd

### DE resulsts and figures for expanded T cell CITE proteins and RNA ###

#TCR data processing
Rscript tcr_processing/A_compile_tcr_data_2020-07-28.R
Rscript tcr_processing/B_define_expansion.R

#pseudobulk_pooling and differential expression CITE
#pooling
bash pseudobulk_pooling/pool_expanded_clone_CITE/sm_call_pooling
#DE 
bash differential_expression/expanded_clone_pseudobulk_CITE_DE/sm_call_t0_covid
bash differential_expression/expanded_clone_pseudobulk_CITE_DE/sm_call_t0_whealthy
bash differential_expression/expanded_clone_pseudobulk_CITE_DE/sm_call_covid_timecourse

#figures CITE
knit figures/tcell_exhaustion/exhaustion_fgsea_plots.Rmd
knit figures/tcell_exhaustion/cutoff_markers_grid.Rmd
knit figures/tcell_exhaustion/CITE_surface_marker_coefficients_with_pval.Rmd

#pseudobulk_pooling and differential expression RNA

bash pseudobulk_pooling/pool_expanded_clone_RNA/sm_call_pooling

bash differential_expression/expanded_clone_pseudobulk_RNA/sm_call_t0_covid
bash differential_expression/expanded_clone_pseudobulk_RNA/sm_call_t0_whealthy
bash differential_expression/expanded_clone_pseudobulk_RNA/sm_call_covid_timecourse

#figures RNA
knit figures/tcell_exhaustion/RNA_surface_marker_coeff_with_pval.Rmd


### cell frequency differential abundance models ###
# run models
bash differential_expression/covid_de_pipeline_cell_freq/sm_call_master
#figures
knit figures/cd4_8_cluster_volcano.Rmd


### Validation with Sculte-Schrepping et al. 2020 ###

#cohort1
#update metadata to add new columns needed for models
Rscript validation_schulteschrepping/update_seurat_meta_cohort1.R
#pooling
bash validation_schulteschrepping/pseudobulk_pooling/cohort1/sm_call_pooling
#DE
bash validation_schulteschrepping/differential_expression/cohort1/sm_call_covid_timecourse
bash validation_schulteschrepping/differential_expression/cohort1/sm_call_t0_covid
bash validation_schulteschrepping/differential_expression/cohort1/sm_call_t0_whealthy

#cohort2
#update metadata to add new columns needed for models
Rscript validation_schulteschrepping/update_seurat_meta_cohort2.R
#pooling
bash validation_schulteschrepping/pseudobulk_pooling/cohort2/sm_call_pooling
#DE
bash validation_schulteschrepping/differential_expression/cohort2/sm_call_covid_timecourse
bash validation_schulteschrepping/differential_expression/cohort2/sm_call_t0_covid
bash validation_schulteschrepping/differential_expression/cohort2/sm_call_t0_whealthy

#figures
knit figures/validations/NK_fattyacid_gsea_curves_cohorts_1_2.Rmd
knit figures/validations/pdc_apoptosis/pdc_apoptosis_schulteschrepping.Rmd
knit figures/validations/exhaustion_plots.Rmd


### T cell clonality ###

#TCR data processing
Rscript tcr_processing/A_compile_tcr_data_2020-07-28.R
# run model
bash tcell_clonality/sm_call
#figures
knit figures/clonality/plot_pcgroup_onset_assoc_for_paper.Rmd

