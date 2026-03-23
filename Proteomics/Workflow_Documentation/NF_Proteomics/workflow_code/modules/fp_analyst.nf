process FP_ANALYST {
    tag "${data_type}"
    containerOptions = "--cleanenv --bind \$PWD,\$HOME/.config,${projectDir}"
    
    publishDir "${output_dir}/FragPipeAnalystR/${data_type}/",
        mode: params.publish_dir_mode,
        pattern: "output/**",
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    publishDir "${output_dir}/FragPipeAnalystR/${data_type}/",
        mode: params.publish_dir_mode,
        pattern: "versions.yml"

    input:
    val(output_dir)
    tuple val(data_type), path(quantification_file), path(experiment_annotation), val(gene_annotations_url)
    // data_type: level for R script (protein, gene, peptide, site). LFQ: protein. TMT: protein/gene/peptide/site.
    // Feature lists: protein/gene for protein level; peptide for peptide level; site for site level. Empty = use top_n_*.

    output:
    path("output/**"), emit: output_files
    path("versions.yml"), emit: versions

    script:
    // Derive mode from fragpipe_workflow parameter (LFQ/TMT)
    def mode = params.fragpipe_workflow?.startsWith('TMT') ? 'TMT' : 'LFQ'
    
    def feature_list_protein = (params.fp_analyst_protein_feature_list != null && params.fp_analyst_protein_feature_list != '') ? params.fp_analyst_protein_feature_list : ''
    def feature_list_gene = (params.fp_analyst_gene_feature_list != null && params.fp_analyst_gene_feature_list != '') ? params.fp_analyst_gene_feature_list : ''
    def feature_list_peptide = (params.fp_analyst_peptide_feature_list != null && params.fp_analyst_peptide_feature_list != '') ? params.fp_analyst_peptide_feature_list : ''
    def feature_list_site = (params.fp_analyst_site_feature_list != null && params.fp_analyst_site_feature_list != '') ? params.fp_analyst_site_feature_list : ''
    def top_n_protein = params.fp_analyst_top_n_protein != null ? params.fp_analyst_top_n_protein : 10
    def top_n_gene = params.fp_analyst_top_n_gene != null ? params.fp_analyst_top_n_gene : 10
    def top_n_peptide = params.fp_analyst_top_n_peptide != null ? params.fp_analyst_top_n_peptide : 10
    def top_n_site = params.fp_analyst_top_n_site != null ? params.fp_analyst_top_n_site : 10
    def enrichment_database = (params.fp_analyst_enrichment_database != null && params.fp_analyst_enrichment_database != '') ?
        (params.fp_analyst_enrichment_database instanceof List ? params.fp_analyst_enrichment_database.join(',') : params.fp_analyst_enrichment_database.toString()) : ''
    def enrichment_direction = (params.fp_analyst_enrichment_direction != null && params.fp_analyst_enrichment_direction != '') ? params.fp_analyst_enrichment_direction : 'Up,Down'
    def gsea_database = (params.fp_analyst_gsea_database != null && params.fp_analyst_gsea_database != '') ?
        (params.fp_analyst_gsea_database instanceof List ? params.fp_analyst_gsea_database.join(',') : params.fp_analyst_gsea_database.toString()) : ''
    def lfq_type = params.fp_analyst_lfq_type ?: "Intensity"
    def normalization_method = (params.fp_analyst_normalization_method != null && params.fp_analyst_normalization_method != '') ? params.fp_analyst_normalization_method : 'none'
    def de_alpha = params.fp_analyst_de_alpha ?: 0.05
    def de_lfc = params.fp_analyst_de_lfc ?: 1.0
    def de_fdr = (params.fp_analyst_de_fdr != null && params.fp_analyst_de_fdr != '') ? params.fp_analyst_de_fdr : 'Benjamini Hochberg'
    def imputation_type = (params.fp_analyst_imputation_type != null && params.fp_analyst_imputation_type != '') ? params.fp_analyst_imputation_type : 'Perseus-type'
    def imputation_shift = params.fp_analyst_imputation_shift != null ? params.fp_analyst_imputation_shift : 1.8
    def imputation_scale = params.fp_analyst_imputation_scale != null ? params.fp_analyst_imputation_scale : 0.3
    // def min_global = params.fp_analyst_min_global_appearance != null ? params.fp_analyst_min_global_appearance : 0
    // def min_cond = params.fp_analyst_min_appearance_one_cond != null ? params.fp_analyst_min_appearance_one_cond : 0
    def qc_plot_data = (params.fp_analyst_qc_plot_data != null && params.fp_analyst_qc_plot_data != '') ? params.fp_analyst_qc_plot_data.toString().toLowerCase() : 'nonimputed'
    def sample_cvs_full_range = (params.fp_analyst_sample_cvs_full_range == true || params.fp_analyst_sample_cvs_full_range == 'true') ? 'true' : 'false'
    def volcano_display_names = (params.fp_analyst_volcano_display_names == true || params.fp_analyst_volcano_display_names == 'true') ? 'true' : 'false'
    def volcano_show_gene = (params.fp_analyst_volcano_show_gene == true || params.fp_analyst_volcano_show_gene == 'true') ? 'true' : 'false'
    def assay_suffix = (params.assay_suffix != null && params.assay_suffix != '') ? params.assay_suffix.toString() : ''
    def gene_annotations_arg = (gene_annotations_url == null || gene_annotations_url?.toString()?.trim() == '') ? '' : "--gene_annotations \"${gene_annotations_url}\""
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
        --enrichment_database "${enrichment_database}" \\
        --enrichment_direction "${enrichment_direction}" \\
        --gsea_database "${gsea_database}" \\
        --lfq_type "${lfq_type}" \\
        --normalization_method "${normalization_method}" \\
        --de_alpha "${de_alpha}" \\
        --de_lfc "${de_lfc}" \\
        --de_fdr "${de_fdr}" \\
        --imputation_type "${imputation_type}" \\
        --imputation_shift "${imputation_shift}" \\
        --imputation_scale "${imputation_scale}" \\
        --qc_plot_data "${qc_plot_data}" \\
        --sample_cvs_full_range "${sample_cvs_full_range}" \\
        --volcano_display_names "${volcano_display_names}" \\
        --volcano_show_gene "${volcano_show_gene}" \\
        ${gene_annotations_arg} \\
        --output_dir "output/"

    # Version info
    Rscript -e "
    versions <- c();
    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
    tryCatch({ versions['FragPipeAnalystR'] <- as.character(packageVersion('FragPipeAnalystR')) }, error = function(e) { versions['FragPipeAnalystR'] <<- 'not installed' });
    cat('"FP_ANALYST":\\n', paste0('    ', names(versions), ': ', versions, collapse='\\n'), '\\n', sep='', file='versions.yml')
    "
    """
}
