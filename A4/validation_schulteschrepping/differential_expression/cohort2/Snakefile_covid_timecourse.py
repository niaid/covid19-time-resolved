workdir: "data/externalData/SchulteSchrepping/differential_expression_cohort2/sample_groups/all_timepoints_covid_only"

FORMULAS = {
        #"days_onset": "~ days_after_onset + age_numeric",
        #"Severity_onset_group_interaction" : "~ 0 + severity_onsetgroup2 + age_numeric",
        "Severity_onsetcontinous_interaction" : "~ group_per_sample:days_after_onset + group_per_sample + days_after_onset + age_numeric"
        }

COEFFICIENTS = {
        #"days_onset": {
        #    "days_after_onset" : "days_after_onset",
        #    "Age" : "age_numeric"},
        #"Severity_onset_group_interaction" : {
        #    "severe-mild_X_juncture-prejuncture" : "(severity_onsetgroup2severe_juncture17_23 - severity_onsetgroup2mild_juncture17_23) - (severity_onsetgroup2severe_prejuncture7_17 - severity_onsetgroup2mild_prejuncture7_17)",
        #    #"juncture-early_in_PC1high": "PC1_onset_groupPC1_high_mid - PC1_onset_groupPC1_high_early",
        #    #"mid-early_in_PC1low": "PC1_onset_groupPC1_low_mid - PC1_onset_groupPC1_low_early",
        #    "severe-mild_juncture" : "severity_onsetgroup2severe_juncture17_23 - severity_onsetgroup2mild_juncture17_23",
        #    "severe-mild_prejuncture" : "severity_onsetgroup2severe_prejuncture7_17 - severity_onsetgroup2mild_prejuncture7_17",
        #    },
        "Severity_onsetcontinous_interaction" : {
            "days_onset_in_severe" : "group_per_samplesevereXdays_after_onset + days_after_onset",
            "days_onset_in_mild" : "days_after_onset",
            "diff_days_onset_severe-mild" : "group_per_samplesevereXdays_after_onset"
            }
            #"Age" : "age_numeric"},
        }

CELLTYPES = []
#dge_files = os.listdir("pseudobulk_dgelists")
import os
from glob import glob
#Recursively listing files
# see https://stackoverflow.com/questions/18394147/recursive-sub-folder-search-and-return-files-in-a-list-python
dge_files = [y for x in os.walk("pseudobulk_dgelists") for y in glob(os.path.join(x[0], '*.rds'))]

for fp in dge_files:
    cell = fp
    cell = cell.replace("-pseudobulk_dgelist.rds", "")
    cell = cell.replace("pseudobulk_dgelists/", "")
    CELLTYPES.append(cell)



#CELLTYPES = CELLTYPES[6:7]

rule all:
    input:
        expand("results/{formula}/design_matrices/{celltype}-model@{formula}-design.rds",
               formula = FORMULAS.keys(), celltype = CELLTYPES),
        expand("results/{formula}/limma_voom_fit/{celltype}-model@{formula}-fit.rds", 
               formula = FORMULAS.keys(), celltype = CELLTYPES),
        [expand("results/{formula}/volcano_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--volcano.pdf", 
               formula = FORM, celltype = CELLTYPES, coef = COEFFICIENTS[FORM].keys()) for FORM in FORMULAS.keys()],
        [expand("results/{formula}/fgsea_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.pdf", 
               formula = FORM, celltype = CELLTYPES, coef = COEFFICIENTS[FORM].keys()) for FORM in FORMULAS.keys()],
        #[expand("results/{formula}/fgsea_le_overlap_by_celltype/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea_overlap.pdf",
        #       formula = FORM, celltype = CELLTYPES, coef = COEFFICIENTS[FORM].keys()) for FORM in FORMULAS.keys()],
        #[expand("results/{formula}/camera_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--camera.pdf", 
               #formula = FORM, celltype = CELLTYPES, coef = COEFFICIENTS[FORM]) for FORM in FORMULAS.keys()]
        
        
#compile_patient_metadata.R
#add_metadata_seurat.R
#make_dgelist_all_cells

rule normalize_dgelist:
    input:
        #"data/CITE5p/all_batches/sorted/time_DE/pseudobulk_dgelists/{celltype}-pseudobulk_dgelist.rds"
        "pseudobulk_dgelists/{celltype}-pseudobulk_dgelist.rds"
    output:
        "pseudobulk_dgelists_normalized/{celltype}-pseudobulk_dgelist_normalized.rds"
    script:
        "limma_scripts/calcnormfactors_and_filter_genes_samples.R"

rule make_design_matrices:
    input:
        "pseudobulk_dgelists_normalized/{celltype}-pseudobulk_dgelist_normalized.rds"
    output:
        "results/{formula}/design_matrices/{celltype}-model@{formula}-design.rds"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.formula]
    script:
        "limma_scripts/make_design_matrix.R"

rule run_de:
    input:
        dgelist="pseudobulk_dgelists_normalized/{celltype}-pseudobulk_dgelist_normalized.rds",
        design="results/{formula}/design_matrices/{celltype}-model@{formula}-design.rds"
    output:
        "results/{formula}/limma_voom_fit/{celltype}-model@{formula}-fit.rds"
    params:
        patient_id_col="Subject"
    script:
        "limma_scripts/run_de.R"

rule make_contrast_matrix:
    input:
        "results/{formula}/design_matrices/{celltype}-model@{formula}-design.rds"
    output:
        "results/{formula}/contrast_matrices/{celltype}--model@{formula}--coef@{coef}--contrast_mat.rds"
    params:
        contrast = lambda wildcards : COEFFICIENTS[wildcards.formula][wildcards.coef]
    script:
        "limma_scripts/make_contrast_matrix.R"

rule contrast_fit:
    input:
        "results/{formula}/limma_voom_fit/{celltype}-model@{formula}-fit.rds",
        "results/{formula}/contrast_matrices/{celltype}--model@{formula}--coef@{coef}--contrast_mat.rds"
    output:
        "results/{formula}/contrast_fit/{celltype}--model@{formula}--coef@{coef}--cfit.rds"
    script:
        "limma_scripts/contrast_fit.R"

rule get_toptab:
    input:
        "results/{formula}/contrast_fit/{celltype}--model@{formula}--coef@{coef}--cfit.rds"
    output:
        "results/{formula}/toptables/{coef}/{celltype}--model@{formula}--coef@{coef}--toptab.tsv"
    script:
        "limma_scripts/get_toptable.R"

rule plot_volcano:
    input:
        "results/{formula}/toptables/{coef}/{celltype}--model@{formula}--coef@{coef}--toptab.tsv"
    output:
        "results/{formula}/volcano_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--volcano.pdf"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.formula]
    script:
        "limma_scripts/plot_volcano.R"

rule camera:
    input:
        "results/{formula}/toptables/{coef}/{celltype}--model@{formula}--coef@{coef}--toptab.tsv",
        "genesets/kegg_go_btm_reactome_foointerferon.rds"
    output:
        "results/{formula}/camera_tables/{coef}/{celltype}--model@{formula}--coef@{coef}--camera.tsv"
    script:
        "run_camera_enrich_snakemake.R"

rule plot_camera:
    input:
        "results/{formula}/camera_tables/{coef}/{celltype}--model@{formula}--coef@{coef}--camera.tsv"
    output:
        "results/{formula}/camera_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--camera.pdf"
    script:
        "plot_enrichment_snakemake.R"

rule fgsea:
    input:
        "results/{formula}/toptables/{coef}/{celltype}--model@{formula}--coef@{coef}--toptab.tsv",
        "genesets/kegg_go_btm_reactome_foointerferon.rds"
    output:
        "results/{formula}/fgsea_tables/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.tsv",
        "results/{formula}/fgsea_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.pdf"
    script:
        "run_fgsea.R"

rule fgsea_le_overlap_by_celltype:
    input:
        "results/{formula}/fgsea_tables/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.tsv",
        "genesets/pathway_summary-reduced.AJM JT.xlsx"
    output:
        "results/{formula}/fgsea_le_overlap_by_celltype/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea_overlap.pdf",
    script:
        "jaccard_le_heatmaps_individual_celltype.R"

