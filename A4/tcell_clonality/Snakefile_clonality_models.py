workdir: "data/CITE5p/all_batches/tcr/diversity_pipeline/2020_07_30/"

KEEP_CELLTYPES = {
        "CD4_Mem" : "CD4_Mem", 
        "CD8_Mem" : "CD8_Mem"
}

MIN_CELLS = {
        "CD4_Mem" : "50", 
        "CD8_Mem" : "50"
}

KEEP_SORTED = {
        "Sorted" : "Y", 
        "Unsorted" : "N"
}

KEEP_BatchES = {
        "Unsorted" : ["B1", "B2", "B3"], 
        "Sorted" : ["B2", "B3"]
}

KEEP_TIMEPOINTS = { 
        "t0_covid_only" : ["T0"],
        "t0_plus_healthy" : ["T0", "HC"],
        "all_timepoints_covid_only" : ["T0", "T1", "T2", "T3"]
        }

SAMPLE_GROUPS = KEEP_TIMEPOINTS.keys()

FORMULAS_T0_WHEALTHY = {
        "healthy_vs_covid": "~ Class + Age  + (1|Batch)",
        }

COEFFICIENTS_T0_WHEALTHY = {
        "healthy_vs_covid" : "ClassCOVID",
        }


FORMULAS_COVID_TIMECOURSE = {
        "days_onset": "~ days_since_onset + Age + (1|Batch)",
        }

COEFFICIENTS_COVID_TIMECOURSE = {
        "days_onset": "days_since_onset",
        }


FORMULAS_T0_COVID = {
        "PC1": "~ PC1 + days_since_onset + Age + (1|Batch)",
        "PC1_cat": "~ PC1_cat + days_since_onset + Age + (1|Batch)",
        }

COEFFICIENTS_T0_COVID = {
        "PC1" : "PC1",
        "PC1_cat" : "PC1_catPC1_high",
        }

COEFFICIENTS = {
        "t0_covid_only" : COEFFICIENTS_T0_COVID,
        "t0_plus_healthy" : COEFFICIENTS_T0_WHEALTHY,
        "all_timepoints_covid_only" : COEFFICIENTS_COVID_TIMECOURSE
        }

FORMULAS = {
        "t0_covid_only" : FORMULAS_T0_COVID,
        "t0_plus_healthy" : FORMULAS_T0_WHEALTHY,
        "all_timepoints_covid_only" : FORMULAS_COVID_TIMECOURSE
        }

rule all:
    input:
        #expand("subsetted_clones/{celltypes}-{sorted}_clones.rds", celltypes = KEEP_CELLTYPES.keys(), sorted = KEEP_SORTED.keys()),
        #expand("cells_per_sample_plots/{celltypes}-{sorted}_cells_per_sample.pdf", celltypes = KEEP_CELLTYPES.keys(), sorted = KEEP_SORTED.keys()),
        #expand("diversity_metrics/{celltypes}-{sorted}_diversity.rds", celltypes = KEEP_CELLTYPES.keys(), sorted = KEEP_SORTED.keys()),
        #"sample_meta.tsv"
        #[[expand("sample_groups/{sample_group}/results/results/{formula}/volcano_plots/{coef}/{cell_annotation}-model@results/{formula}-coef@{coef}-volcano.pdf", 
        #       formula = FORM, cell_annotation = CELL_ANNOTATION, coef = COEFFICIENTS[SAMPG][FORM], sample_group = SAMPG) for FORM in FORMULAS[SAMPG].keys()] for SAMPG in SAMPLE_GROUPS]
        [[expand("sample_groups/{sample_group}/results/{formula}/lmertest_summary/{celltypes}-{sorted}_summary_dat.tsv",
               formula = FORM, celltypes = KEEP_CELLTYPES.keys(), sorted = KEEP_SORTED.keys(), sample_group = SAMPG) for FORM in FORMULAS[SAMPG].keys()] for SAMPG in SAMPLE_GROUPS]

        
rule subset_clones:
    #Subset the clonotype data by celltype and sorting 
    # remove Batch 1 from the sorted data
    input:
        "data/CITE5p/all_batches/tcr/2020_07_28_tenx_filtered_anno_wfiltered_celltype.rds"
    output:
        "subsetted_clones/{celltypes}-{sorted}_clones.rds"
    params:
        keep_celltypes = lambda wildcards: KEEP_CELLTYPES[wildcards.celltypes],
        keep_Batches = lambda wildcards: KEEP_BatchES[wildcards.sorted],
        keep_sorted = lambda wildcards: KEEP_SORTED[wildcards.sorted]
    script:
        "subset_clones.R"
   
rule plot_cells_per_sample:
    input:
        "subsetted_clones/{celltypes}-{sorted}_clones.rds"
    output:
        "cells_per_sample_plots/{celltypes}-{sorted}_cells_per_sample.pdf"
    params:
        n_cell_cutoff = 100
    script:
        "plot_cells_per_sample.R"

rule calc_diversity:
    input:
        "subsetted_clones/{celltypes}-{sorted}_clones.rds"
    output:
        "diversity_metrics/{celltypes}-{sorted}_diversity.rds"
    params:
        n_cell_cutoff = 100
    script:
        "calc_diversity.R"

#rule plot_diversity:
#    input:
#        "diversity_metrics/{celltypes}-{sorted}_diversity.rds"
#    output:
#        "diversity_plots/{celltypes}-{sorted}_diversity.pdf"

rule get_sample_meta:
    #Collapse cell level metadata with updated features to sample level meta with desired features for modeling
    input:
        #This is from update_seurat_meta.R in covid_de_pipeline
        "data/CITE5p/all_batches/2020_08_09.rmBuffSHD8.allcelltypelabels.merge.SNG.wmeta.WithinBatchClustered_metadata.Rds"
    output:
        "sample_meta.tsv"
    script:
        "get_sample_meta.R"

rule subset_diversity:
    #Subset the diversity metrics DF to each sample group
    input:
        "diversity_metrics/{celltypes}-{sorted}_diversity.rds"
    output:
        "sample_groups/{sample_group}/diversity_metrics/{celltypes}-{sorted}_diversity.rds"
    params:
        keep_timepoints = lambda wildcards: KEEP_TIMEPOINTS[wildcards.sample_group]
    script:
        "subset_diversity.R"

rule make_dat_list:
    # Add the metadata to the diversity metric dataframes and put into list. One DF per diversity metric
    input:
        "sample_groups/{sample_group}/diversity_metrics/{celltypes}-{sorted}_diversity.rds",
        "sample_meta.tsv"
    output:
        "sample_groups/{sample_group}/dat_list/{celltypes}-{sorted}_diversity.rds"
    script:
        "make_dat_list.R"

rule filter_dat_list:
    #remove the samples with NA values for any of the things in the formula
    input:
        "sample_groups/{sample_group}/dat_list/{celltypes}-{sorted}_diversity.rds"
    output:
        "sample_groups/{sample_group}/results/{formula}/filtered_dat_list/{celltypes}-{sorted}_diversity.rds"
    params:
        formula = lambda wildcards: FORMULAS[wildcards.sample_group][wildcards.formula]
    script:
        "filter_dat_list.R"

rule run_models:
    input:
        "sample_groups/{sample_group}/results/{formula}/filtered_dat_list/{celltypes}-{sorted}_diversity.rds"
    output:
        "sample_groups/{sample_group}/results/{formula}/lmer_fit/{celltypes}-{sorted}_fit_list.rds"
    params:
        formula = lambda wildcards: FORMULAS[wildcards.sample_group][wildcards.formula],
        y_col = "median1000"
    script:
        "fit_model.R"

rule do_summary:
    input:
        "sample_groups/{sample_group}/results/{formula}/lmer_fit/{celltypes}-{sorted}_fit_list.rds"
    output:
        "sample_groups/{sample_group}/results/{formula}/lmertest_summary/{celltypes}-{sorted}_summary_dat.tsv"
    script:
        "get_summary_dat.R"

