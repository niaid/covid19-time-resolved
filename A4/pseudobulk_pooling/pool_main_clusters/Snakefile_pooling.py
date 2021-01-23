workdir: "data/CITE5p/all_batches/differential_expression/2020_08_09/"

SORT_GROUPS = ["Sorted", "Unsorted"]

CELL_ANNOTATION = ["specific_gating", "WCTcoursecelltype"]

LIBSIZE_FILTER = 40000
NCELL_FILTER = 7


KEEP_TIMEPOINTS = { 
        "t0_covid_only" : ["T0"],
        "t0_plus_healthy" : ["T0", "HC"],
        "all_timepoints_covid_only" : ["T0", "T1", "T2", "T3"],
        #"all_samples" : ["HC", "T0", "T1", "T2", "T3"]
        }
SAMPLE_GROUPS = KEEP_TIMEPOINTS.keys()

KEEP_SORTED = {
        "Sorted" : "Y",
        "Unsorted" : "N" 
        }

KEEP_BATCHES = {
        "Sorted" : ["B2", "B3"],
        "Unsorted" : ["B1", "B2", "B3"]
        }

MIN_SAMPLES_PER_CELLTYPE = {"gate_group_merged" : 5, 
        "adjustedcelltype" : 5, 
        "specific_gating" : 5,
        "WCTcoursecelltype" : 5,
        "WCTmergedcelltype" : 5}

rule all:
    input:
        expand("sample_groups/{sample_group}/pseudobulk_dgelists/{sort_group}-{cell_annotation}/", sample_group = SAMPLE_GROUPS, sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION),
        expand("sample_groups/{sample_group}/pseudobulk_stats_prefilter/{sort_group}-{cell_annotation}.pdf", sample_group = SAMPLE_GROUPS, sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION),
        expand("sample_groups/{sample_group}/pseudobulk_stats_postfilter/{sort_group}-{cell_annotation}.pdf", sample_group = SAMPLE_GROUPS, sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION),
        expand("all_samples/pseudobulk_stats_prefilter/{sort_group}-{cell_annotation}.pdf", sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION),
        expand("all_samples/severity_network_pseudobulk_dgelists_filtered_samples_genes/{sort_group}-{cell_annotation}.rds", sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION),
        expand("all_samples/severity_network_pseudobulk_dgelists_filtered_samples/{sort_group}-{cell_annotation}.rds", sort_group = SORT_GROUPS, cell_annotation = CELL_ANNOTATION)
        
        
rule pseudobulk_pool:
    input:
        seurat="input/brescia_paper1_seurat.rds"
    output:
        "all_samples/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    params:
        keep_batches = lambda wildcards : KEEP_BATCHES[wildcards.sort_group],
        keep_sorted = lambda wildcards : KEEP_SORTED[wildcards.sort_group],
        cell_anno_col = lambda wildcards : wildcards.cell_annotation,
        libsize_filter=1,
        min_cells_per_pool=1,
        min_samples_per_celltype=1
    script:
        "pseudobulk_pooling_snakemake.R"

rule subset_pseudobulk:
    input:
        "all_samples/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        "sample_groups/{sample_group}/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    params:
        keep_timepoints = lambda wildcards : KEEP_TIMEPOINTS[wildcards.sample_group]
    script:
        "subset_pseudobulk_timepoint.R"


rule filter_pseudobulk_list:
    input:
        "sample_groups/{sample_group}/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        directory("sample_groups/{sample_group}/pseudobulk_dgelists/{sort_group}-{cell_annotation}/")
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
        min_samples_per_celltype=lambda wildcards : MIN_SAMPLES_PER_CELLTYPE[wildcards.cell_annotation]
    script:
        "filter_pseudobulk_list.R"

rule plot_pseudobulk_stats_prefilter_all_samples:
    input:
        "all_samples/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        "all_samples/pseudobulk_stats_prefilter/{sort_group}-{cell_annotation}.pdf"
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
    script:
        "pseudobulk_stats.R"

rule plot_pseudobulk_stats_prefilter:
    input:
        "sample_groups/{sample_group}/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        "sample_groups/{sample_group}/pseudobulk_stats_prefilter/{sort_group}-{cell_annotation}.pdf"
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
    script:
        "pseudobulk_stats.R"

rule plot_pseudobulk_stats_post_filter:
    input:
        "sample_groups/{sample_group}/pseudobulk_dgelists/{sort_group}-{cell_annotation}/"
    output:
        "sample_groups/{sample_group}/pseudobulk_stats_postfilter/{sort_group}-{cell_annotation}.pdf"
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
    script:
        "pseudobulk_stats.R"

rule filter_pseudobulk_for_severity_network_samples_only:
    input:
        "all_samples/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        "all_samples/severity_network_pseudobulk_dgelists_filtered_samples/{sort_group}-{cell_annotation}.rds"
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
        min_samples_per_celltype=lambda wildcards : MIN_SAMPLES_PER_CELLTYPE[wildcards.cell_annotation]
    script:
        "filter_pseudobulk_list.R"

rule filter_pseudobulk_for_severity_network:
    input:
        "all_samples/pseudobulk_dgelists_unfiltered/{sort_group}-{cell_annotation}.rds"
    output:
        "all_samples/severity_network_pseudobulk_dgelists_filtered_samples_genes/{sort_group}-{cell_annotation}.rds"
    params:
        libsize_filter=LIBSIZE_FILTER,
        min_cells_per_pool=NCELL_FILTER,
        min_samples_per_celltype=lambda wildcards : MIN_SAMPLES_PER_CELLTYPE[wildcards.cell_annotation]
    script:
        "filter_pseudobulk_for_severity_network.R"

