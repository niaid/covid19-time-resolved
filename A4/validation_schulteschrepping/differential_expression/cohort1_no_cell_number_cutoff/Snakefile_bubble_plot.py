rule fgsea_bubble:
    input:
        "results/{formula}/fgsea_tables/{coef}/{sorted}-{cell_anno}/",
    output:
        "results/{formula}/fgsea_bubble_plots/{coef}/{sorted}-{cell_anno}_btm-kegg.png",
        "results/{formula}/fgsea_bubble_plots/{coef}/{sorted}-{cell_anno}_go.png",
        "results/{formula}/fgsea_bubble_plots/{coef}/{sorted}-{cell_anno}_reactome.png"
    script:
        "plotting/plot_fgsea_bubble.R"
