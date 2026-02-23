process FP_ANALYST {
    tag "${data_type}"
    containerOptions = "--cleanenv --bind \$PWD,\$HOME/.config,${projectDir}"
    
    publishDir "${output_dir}/FragPipe-Analyst/${data_type}/",
        mode: params.publish_dir_mode,
        pattern: "output/**",
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    publishDir "${output_dir}/FragPipe-Analyst/${data_type}/",
        mode: params.publish_dir_mode,
        pattern: "versions.yml"

    input:
    val(output_dir)
    tuple val(data_type), path(quantification_file), path(experiment_annotation)
    // data_type: level for R script (protein, gene, peptide, site). LFQ: protein. TMT: protein/gene/peptide/site.
    // Feature lists: protein/gene for protein level; peptide for peptide level; site for site level. Empty = use top_n_*.
    // Pathway/GO params: fp_analyst_pathway_database, fp_analyst_pathway_direction, fp_analyst_go_database, fp_analyst_go_direction

    output:
    path("output/**"), emit: output_files
    path("versions.yml"), emit: versions

    script:
    // Derive mode from fragpipe_workflow parameter
    def mode = params.fragpipe_workflow?.startsWith('TMT') ? 'TMT' : 
               params.fragpipe_workflow?.startsWith('DIA') ? 'DIA' : 
               'LFQ'
    
    def feature_list_protein = (params.fp_analyst_protein_feature_list != null && params.fp_analyst_protein_feature_list != '') ? params.fp_analyst_protein_feature_list : ''
    def feature_list_gene = (params.fp_analyst_gene_feature_list != null && params.fp_analyst_gene_feature_list != '') ? params.fp_analyst_gene_feature_list : ''
    def feature_list_peptide = (params.fp_analyst_peptide_feature_list != null && params.fp_analyst_peptide_feature_list != '') ? params.fp_analyst_peptide_feature_list : ''
    def feature_list_site = (params.fp_analyst_site_feature_list != null && params.fp_analyst_site_feature_list != '') ? params.fp_analyst_site_feature_list : ''
    def top_n_protein = params.fp_analyst_top_n_protein != null ? params.fp_analyst_top_n_protein : 10
    def top_n_gene = params.fp_analyst_top_n_gene != null ? params.fp_analyst_top_n_gene : 10
    def top_n_peptide = params.fp_analyst_top_n_peptide != null ? params.fp_analyst_top_n_peptide : 10
    def top_n_site = params.fp_analyst_top_n_site != null ? params.fp_analyst_top_n_site : 10
    def pathway_database = (params.fp_analyst_pathway_database != null && params.fp_analyst_pathway_database != '') ?
        (params.fp_analyst_pathway_database instanceof List ? params.fp_analyst_pathway_database.join(',') : params.fp_analyst_pathway_database.toString()) : ''
    def pathway_direction = (params.fp_analyst_pathway_direction != null && params.fp_analyst_pathway_direction != '') ? params.fp_analyst_pathway_direction : 'Both'
    def go_database = (params.fp_analyst_go_database != null && params.fp_analyst_go_database != '') ?
        (params.fp_analyst_go_database instanceof List ? params.fp_analyst_go_database.join(',') : params.fp_analyst_go_database.toString()) : ''
    def go_direction = (params.fp_analyst_go_direction != null && params.fp_analyst_go_direction != '') ? params.fp_analyst_go_direction : 'Both'
    def lfq_type = params.fp_analyst_lfq_type ?: "Intensity"
    def normalization_method = (params.fp_analyst_normalization_method != null && params.fp_analyst_normalization_method != '') ? params.fp_analyst_normalization_method : 'none'
    def de_alpha = params.fp_analyst_de_alpha ?: 0.05
    def de_lfc = params.fp_analyst_de_lfc ?: 1.0
    def de_fdr = (params.fp_analyst_de_fdr != null && params.fp_analyst_de_fdr != '') ? params.fp_analyst_de_fdr : 'Benjamini Hochberg'
    def imputation_type = (params.fp_analyst_imputation_type != null && params.fp_analyst_imputation_type != '') ? params.fp_analyst_imputation_type : 'Perseus-type'
    def min_global = params.fp_analyst_min_global_appearance != null ? params.fp_analyst_min_global_appearance : 0
    def min_cond = params.fp_analyst_min_appearance_one_cond != null ? params.fp_analyst_min_appearance_one_cond : 0
    def qc_imputed = (params.fp_analyst_qc_show_imputed == true || params.fp_analyst_qc_show_imputed == 'true') ? 'true' : 'false'
    def qc_both = (params.fp_analyst_qc_include_both == true || params.fp_analyst_qc_include_both == 'true') ? 'true' : 'false'
    def volcano_display_names = (params.fp_analyst_volcano_display_names == true || params.fp_analyst_volcano_display_names == 'true') ? 'true' : 'false'
    def volcano_show_gene = (params.fp_analyst_volcano_show_gene == true || params.fp_analyst_volcano_show_gene == 'true') ? 'true' : 'false'
    def volcano_highlight_feature = (params.fp_analyst_volcano_highlight_feature != null && params.fp_analyst_volcano_highlight_feature != '') ? params.fp_analyst_volcano_highlight_feature : ''
    def volcano_show_other_peptides = (params.fp_analyst_volcano_show_other_peptides == true || params.fp_analyst_volcano_show_other_peptides == 'true') ? 'true' : 'false'
    def pkg_list = ['SummarizedExperiment', 'dplyr', 'tibble', 'tidyr', 'purrr', 'ggplot2', 'matrixStats', 'limma', 'ComplexHeatmap', 'circlize', 'RColorBrewer', 'ggrepel', 'scales', 'vegan', 'cluster', 'httr', 'data.table', 'ggVennDiagram', 'UpSetR']
    if (imputation_type != 'Perseus-type' && imputation_type != 'none') pkg_list << 'MSnbase'
    if (normalization_method == 'vsn') pkg_list << 'vsn'
    if (de_fdr != null && !de_fdr.toLowerCase().contains('benjamini') && !de_fdr.equalsIgnoreCase('bh')) pkg_list << 'fdrtool'
    def pkg_list_r = "c(" + pkg_list.collect { "'$it'" }.join(", ") + ")"
    """
    # Create output directory
    mkdir -p output/

    # Run FragPipe-Analyst R script (executable, matches test_fp_analyst_datasets.R)
    fp_analyst_main.R \\
        --experiment_annotation "${experiment_annotation}" \\
        --quantification_file "${quantification_file}" \\
        --mode "${mode}" \\
        --level "${data_type}" \\
        --feature_list_protein "${feature_list_protein}" \\
        --feature_list_gene "${feature_list_gene}" \\
        --feature_list_peptide "${feature_list_peptide}" \\
        --feature_list_site "${feature_list_site}" \\
        --top_n_protein "${top_n_protein}" \\
        --top_n_gene "${top_n_gene}" \\
        --top_n_peptide "${top_n_peptide}" \\
        --top_n_site "${top_n_site}" \\
        --pathway_database "${pathway_database}" \\
        --pathway_direction "${pathway_direction}" \\
        --go_database "${go_database}" \\
        --go_direction "${go_direction}" \\
        --lfq_type "${lfq_type}" \\
        --normalization_method "${normalization_method}" \\
        --de_alpha "${de_alpha}" \\
        --de_lfc "${de_lfc}" \\
        --de_fdr "${de_fdr}" \\
        --imputation_type "${imputation_type}" \\
        --min_global_appearance "${min_global}" \\
        --min_appearance_one_condition "${min_cond}" \\
        --qc_show_imputed "${qc_imputed}" \\
        --qc_include_both "${qc_both}" \\
        --volcano_display_names "${volcano_display_names}" \\
        --volcano_show_gene "${volcano_show_gene}" \\
        --volcano_highlight_feature "${volcano_highlight_feature}" \\
        --volcano_show_other_peptides "${volcano_show_other_peptides}" \\
        --output_dir "output/"
    
    # Version info
    Rscript -e "
    versions <- c();
    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
    tryCatch({ versions['BiocManager'] <- as.character(BiocManager::version()) }, error = function(e) { versions['BiocManager'] <<- 'unknown' });
    pkg_list <- ${pkg_list_r};
    for (pkg in pkg_list) {
        tryCatch(
            { versions[pkg] <- as.character(packageVersion(pkg)) },
            error = function(e) { versions[pkg] <<- 'not installed' }
        )
    };
    cat('"FP_ANALYST":\\n', paste0('    ', names(versions), ': ', versions, collapse='\\n'), '\\n', sep='', file='versions.yml')
    "
    """
}
