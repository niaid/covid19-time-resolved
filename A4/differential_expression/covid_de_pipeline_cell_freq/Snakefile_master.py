workdir: "data/CITE5p/all_batches/differential_expression_cell_freq/2020_08_25/"

SAMPLE_GROUPS = ["t0_covid_only",
        "t0_plus_healthy",
        "all_timepoints_covid_only"]

AUTO_CLUSTER_DIR = "/hpcdata/sg/sg_data/users/liuc19/COVID19/CITE5p_Brescia_3batches/auto_clusters/"

MAN_GATE_DIR = "/hpcdata/sg/sg_data/users/liuc19/COVID19/CITE5p_Brescia_3batches/flowjo_gating/"

INPUT_MAT_PATHS = {
        "final_full_UnSort_rmLowParent_20200818" : AUTO_CLUSTER_DIR + "final_full_cell_freq_UnSort_mtx.20200818.csv",
        "TotalCD4and8_LE_genes_20200824" : "data/CITE5p/all_batches/cell_freqs/TotalCD4and8.LG.20200824.merge_cell_UnSort_WCTmerged_to_parent_mtx.Rds"
        }

CELL_ANNOTATION = INPUT_MAT_PATHS.keys()

KEEP_TIMEPOINTS = { 
        "t0_covid_only" : ["T0"],
        "t0_plus_healthy" : ["T0", "HC"],
        "all_timepoints_covid_only" : ["T0", "T1", "T2", "T3"]
        }

FORMULAS_T0_WHEALTHY = {
        "healthy_vs_covid" : "~ 0 + Class + Age + Batch",
        }

COEFFICIENTS_T0_WHEALTHY = {
        "healthy_vs_covid" : {
            "COVID-Healthy" : "ClassCOVID - ClassHC",
            "Age" : "Age"},
        }

FORMULAS_T0_COVID = {
        "days_since_onset": "~ days_since_onset + Age + batch",
        "PC1": "~ PC1 + days_since_onset + Age + batch",
        "PC1_group": "~ 0 + PC1_cat + days_since_onset + Age + batch",
        }

COEFFICIENTS_T0_COVID = {
        "days_since_onset": {
            "days_since_onset" : "days_since_onset",
            "Age" : "Age"},
        "PC1" : {
            "PC1" : "PC1",
            "Age" : "Age"},
        "PC1_group" : {
            "PC1high-low" : "PC1_catPC1_high - PC1_catPC1_low",
            "Age" : "Age"},
        }

FORMULAS_COVID_TIMECOURSE = {
        "days_onset": "~ days_since_onset + Age + batch",
        "onset_group" : "~ 0 + onset_group + Age + batch",
        "PC1_onset_group_interaction" : "~ 0 + PC1_onset_group + Age + batch",
        "PC1group_onsetcontinous_interaction" : "~ PC1_cat:days_since_onset + PC1_cat + days_since_onset + Age + batch"
        }

COEFFICIENTS_COVID_TIMECOURSE = {
        "days_onset": {
            "days_since_onset" : "days_since_onset",
            "Age" : "Age"},
        "onset_group" : {
            "mid-early": "onset_groupmid - onset_groupearly",
            "mid-late": "onset_groupmid - onset_grouplate",
            "late-early" : "onset_grouplate - onset_groupearly",
            "mid-early&late" : "onset_groupmid - (onset_groupearly + onset_grouplate)/2",
            "Age" : "Age"},
        "PC1_onset_group_interaction" : {
            "PC1High-low_X_onsetMid-early" : "(PC1_onset_groupPC1_high_mid - PC1_onset_groupPC1_low_mid) - (PC1_onset_groupPC1_high_early - PC1_onset_groupPC1_low_early)",
            "mid-early_in_PC1high": "PC1_onset_groupPC1_high_mid - PC1_onset_groupPC1_high_early",
            "mid-early_in_PC1low": "PC1_onset_groupPC1_low_mid - PC1_onset_groupPC1_low_early",
            "PC1High-low_in_mid" : "PC1_onset_groupPC1_high_mid - PC1_onset_groupPC1_low_mid",
            "PC1High-low_in_early" : "PC1_onset_groupPC1_high_early - PC1_onset_groupPC1_low_early",
            "Age" : "Age"},
        "PC1group_onsetcontinous_interaction" : {
            "days_onset_in_PC1high" : "PC1_catPC1_highXdays_since_onset + days_since_onset",
            "days_onset_in_PC1low" : "days_since_onset",
            "diff_days_onset_PC1High-low" : "PC1_catPC1_highXdays_since_onset"
            }
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
        [[expand("{sample_group}/results/{formula}/volcano_plots/{coef}/{cell_annotation}--model@{formula}--coef@{coef}--volcano.pdf", 
               formula = FORM, cell_annotation = CELL_ANNOTATION, coef = COEFFICIENTS[SAMPG][FORM], sample_group = SAMPG) for FORM in FORMULAS[SAMPG].keys()] for SAMPG in SAMPLE_GROUPS]
        
rule create_expressionsets:
    input:
        lambda wildcards: INPUT_MAT_PATHS[wildcards.cell_annotation],
        meta="data/CITE5p/all_batches/2020_08_09.rmBuffSHD8.allcelltypelabels.merge.SNG.wmeta.WithinBatchClustered_metadata.Rds"
    output:
        "expressionsets/{cell_annotation}-eset.rds"
    script:
        "create_eset.R"

rule subset_expressionset:
    input:
        "expressionsets/{cell_annotation}-eset.rds"
    output:
        "{sample_group}/subsetted_expressionsets/{cell_annotation}-eset.rds"
    params:
        keep_timepoints = lambda wildcards : KEEP_TIMEPOINTS[wildcards.sample_group]
    script:
        "subset_expressionset.R"

rule make_design_matrices:
    input:
        "{sample_group}/subsetted_expressionsets/{cell_annotation}-eset.rds"
    output:
        "{sample_group}/results/{formula}/design_matrices/{cell_annotation}-model@{formula}-design.rds"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.sample_group][wildcards.formula]
    script:
        "limma_scripts/make_design_matrix.R"

rule run_de:
    input:
        eset="{sample_group}/subsetted_expressionsets/{cell_annotation}-eset.rds",
        design="{sample_group}/results/{formula}/design_matrices/{cell_annotation}-model@{formula}-design.rds"
    output:
        "{sample_group}/results/{formula}/limma_voom_fit/{cell_annotation}-model@{formula}-fit.rds"
    params:
        patient_id_col="Subject"
    script:
        "limma_scripts/run_de.R"

rule make_contrast_matrix:
    input:
        "{sample_group}/results/{formula}/design_matrices/{cell_annotation}-model@{formula}-design.rds"
    output:
        "{sample_group}/results/{formula}/contrast_matrices/{cell_annotation}--model@{formula}--coef@{coef}--contrast_mat.rds"
    params:
        contrast = lambda wildcards : COEFFICIENTS[wildcards.sample_group][wildcards.formula][wildcards.coef]
    script:
        "limma_scripts/make_contrast_matrix.R"

rule contrast_fit:
    input:
        "{sample_group}/results/{formula}/limma_voom_fit/{cell_annotation}-model@{formula}-fit.rds",
        "{sample_group}/results/{formula}/contrast_matrices/{cell_annotation}--model@{formula}--coef@{coef}--contrast_mat.rds"
    output:
        "{sample_group}/results/{formula}/contrast_fit/{cell_annotation}--model@{formula}--coef@{coef}--cfit.rds"
    script:
        "limma_scripts/contrast_fit.R"

rule get_toptab:
    input:
        "{sample_group}/results/{formula}/contrast_fit/{cell_annotation}--model@{formula}--coef@{coef}--cfit.rds"
    output:
        "{sample_group}/results/{formula}/toptables/{coef}/{cell_annotation}--model@{formula}--coef@{coef}--toptab.tsv"
    script:
        "limma_scripts/get_toptable.R"

rule plot_volcano:
    input:
        "{sample_group}/results/{formula}/toptables/{coef}/{cell_annotation}--model@{formula}--coef@{coef}--toptab.tsv"
    output:
        "{sample_group}/results/{formula}/volcano_plots/{coef}/{cell_annotation}--model@{formula}--coef@{coef}--volcano.pdf"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.sample_group][wildcards.formula]
    script:
        "limma_scripts/plot_volcano.R"
