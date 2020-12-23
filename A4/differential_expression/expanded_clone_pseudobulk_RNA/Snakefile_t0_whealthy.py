workdir: "data/CITE5p/all_batches/expanded_tcell_pbulk_DE/2020_08_25/sample_groups/t0_plus_healthy/"

FORMULAS = {
        "healthy_vs_covid" : "~ 0 + Class + Age + Batch",
        "healthy_vs_severitygroup" : "~ 0 + cond_group + Age + Batch",
        }

COEFFICIENTS = {
        "healthy_vs_covid" : {
            "COVID-Healthy" : "ClassCOVID - ClassHC",
            "Age" : "Age"},
        "healthy_vs_severitygroup" : {
            "Severe-Healthy" : "cond_groupSevere - cond_groupHC",
            "Critical-Healthy" : "cond_groupCritical - cond_groupHC"
            }
        }


#import os
#os.chdir("data/CITE5p/all_batches/differential_expression/2020_07_10/t0_plus_healthy/")

CELLTYPES = []
#dge_files = os.listdir("pseudobulk_dgelists")
import os
from glob import glob
#Recursively listing files
# see https://stackoverflow.com/questions/18394147/recursive-sub-folder-search-and-return-files-in-a-list-python
dge_files = [y for x in os.walk("pseudobulk_dgelists") for y in glob(os.path.join(x[0], '*.rds'))]

for fp in dge_files:
    #cell = os.path.basename(fp)
    cell = fp
    cell = cell.replace("-pseudobulk_dgelist.rds", "")
    cell = cell.replace("pseudobulk_dgelists/", "")
    CELLTYPES.append(cell)



#CELLTYPES = CELLTYPES[0:1]

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

rule fgsea:
    input:
        "results/{formula}/toptables/{coef}/{celltype}--model@{formula}--coef@{coef}--toptab.tsv",
        "genesets/kegg_go_btm_reactome_foointerferon.rds"
    output:
        "results/{formula}/fgsea_tables/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.tsv",
        "results/{formula}/fgsea_plots/{coef}/{celltype}--model@{formula}--coef@{coef}--fgsea.pdf"
    script:
        "run_fgsea.R"
