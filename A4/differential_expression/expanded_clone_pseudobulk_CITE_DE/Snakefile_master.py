workdir: "data/CITE5p/all_batches/expanded_tcell_pbulk_CITE_DE/2020_07_31/sample_groups/"

SAMPLE_GROUPS = ["t0_covid_only",
        "t0_plus_healthy",
        "all_timepoints_covid_only"]

KEEP_TIMEPOINTS = { 
        "t0_covid_only" : ["T0"],
        "t0_plus_healthy" : ["T0", "HC"],
        "all_timepoints_covid_only" : ["T0", "T1", "T2", "T3"]
        }


CELLTYPES = []
#dge_files = os.listdir("pseudobulk_dgelists")
import os
from glob import glob
#Recursively listing files
# see https://stackoverflow.com/questions/18394147/recursive-sub-folder-search-and-return-files-in-a-list-python
#os.chdir("data/CITE5p/all_batches/expanded_tcell_pbulk_CITE_DE/2020_07_31/sample_groups/t0_plus_healthy/")
dge_files = [y for x in os.walk("t0_plus_healthy/pseudobulk_esets") for y in glob(os.path.join(x[0], '*.rds'))]

for fp in dge_files:
    cell = fp
    cell = cell.replace("-pseudobulk_eset.rds", "")
    cell = cell.replace("t0_plus_healthy/pseudobulk_esets/", "")
    CELLTYPES.append(cell)

FORMULAS_T0_WHEALTHY = {
        "healthy_vs_covid": "~ Class + Age  + Batch",
        }

COEFFICIENTS_T0_WHEALTHY = {
        "healthy_vs_covid" : "ClassCOVID",
        }

FORMULAS_COVID_TIMECOURSE = {
        "days_onset": "~ days_since_onset + Age + batch",
        }

COEFFICIENTS_COVID_TIMECOURSE = {
        "days_onset": "days_since_onset",
        }


FORMULAS_T0_COVID = {
        "PC1": "~ PC1 + days_since_onset + Age + batch",
        "PC1_cat": "~ PC1_cat + days_since_onset + Age + batch",
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
        [[expand("{sample_group}/results/{formula}/volcano_plots/{coef}/{cell_annotation}-model@{formula}-coef@{coef}-volcano.pdf", 
               formula = FORM, cell_annotation = CELLTYPES, coef = COEFFICIENTS[SAMPG][FORM], sample_group = SAMPG) for FORM in FORMULAS[SAMPG].keys()] for SAMPG in SAMPLE_GROUPS]
        

rule make_design_matrices:
    input:
        "{sample_group}/pseudobulk_esets/{cell_annotation}-pseudobulk_eset.rds"
    output:
        "{sample_group}/results/{formula}/design_matrices/{cell_annotation}-model@{formula}-design.rds"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.sample_group][wildcards.formula]
    script:
        "limma_scripts/make_design_matrix.R"

rule run_de:
    input:
        eset="{sample_group}/pseudobulk_esets/{cell_annotation}-pseudobulk_eset.rds",
        design="{sample_group}/results/{formula}/design_matrices/{cell_annotation}-model@{formula}-design.rds"
    output:
        "{sample_group}/results/{formula}/limma_voom_fit/{cell_annotation}-model@{formula}-fit.rds"
    params:
        patient_id_col="Subject"
    script:
        "limma_scripts/run_de.R"

rule get_toptab:
    input:
        "{sample_group}/results/{formula}/limma_voom_fit/{cell_annotation}-model@{formula}-fit.rds"
    output:
        "{sample_group}/results/{formula}/toptables/{coef}/{cell_annotation}-model@{formula}-coef@{coef}-toptab.tsv"
    script:
        "limma_scripts/get_toptable.R"

rule plot_volcano:
    input:
        "{sample_group}/results/{formula}/toptables/{coef}/{cell_annotation}-model@{formula}-coef@{coef}-toptab.tsv"
    output:
        "{sample_group}/results/{formula}/volcano_plots/{coef}/{cell_annotation}-model@{formula}-coef@{coef}-volcano.pdf"
    params:
        formula = lambda wildcards : FORMULAS[wildcards.sample_group][wildcards.formula]
    script:
        "limma_scripts/plot_volcano.R"
