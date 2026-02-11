process FP_ANALYST_R {
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
    // Feature lists (both used in single run): protein list plots by rownames, gene list plots by Gene column
    // When a list is empty, top_n_protein or top_n_gene is used to plot top N variable features
    // Pathway/GO params: fp_analyst_pathway_database, fp_analyst_pathway_direction, fp_analyst_go_database, fp_analyst_go_direction

    output:
    path("output/**"), emit: output_files
    path("versions.yml"), emit: versions

    script:
    def script_path = "${projectDir}/bin/fp_analyst_analysis.R"
    
    // Derive mode from fragpipe_workflow parameter
    def mode = params.fragpipe_workflow?.startsWith('TMT') ? 'TMT' : 
               params.fragpipe_workflow?.startsWith('DIA') ? 'DIA' : 
               'LFQ'
    
    def feature_list_protein = (params.fp_analyst_protein_feature_list != null && params.fp_analyst_protein_feature_list != '') ? params.fp_analyst_protein_feature_list : ''
    def feature_list_gene = (params.fp_analyst_gene_feature_list != null && params.fp_analyst_gene_feature_list != '') ? params.fp_analyst_gene_feature_list : ''
    def top_n_protein = params.fp_analyst_top_n_protein != null ? params.fp_analyst_top_n_protein : 10
    def top_n_gene = params.fp_analyst_top_n_gene != null ? params.fp_analyst_top_n_gene : 10
    def pathway_database = (params.fp_analyst_pathway_database != null && params.fp_analyst_pathway_database != '') ? params.fp_analyst_pathway_database : ''
    def pathway_direction = (params.fp_analyst_pathway_direction != null && params.fp_analyst_pathway_direction != '') ? params.fp_analyst_pathway_direction : 'Both'
    def go_database = (params.fp_analyst_go_database != null && params.fp_analyst_go_database != '') ? params.fp_analyst_go_database : ''
    def go_direction = (params.fp_analyst_go_direction != null && params.fp_analyst_go_direction != '') ? params.fp_analyst_go_direction : 'Both'
    def lfq_type = params.fp_analyst_lfq_type ?: "Intensity"
    def de_alpha = params.fp_analyst_de_alpha ?: 0.05
    def de_lfc = params.fp_analyst_de_lfc ?: 1.0
    def de_fdr = (params.fp_analyst_de_fdr != null && params.fp_analyst_de_fdr != '') ? params.fp_analyst_de_fdr : 'Benjamini Hochberg'
    def min_global = params.fp_analyst_min_global_appearance != null ? params.fp_analyst_min_global_appearance : 0
    def min_cond = params.fp_analyst_min_appearance_one_cond != null ? params.fp_analyst_min_appearance_one_cond : 0
    def qc_imputed = (params.fp_analyst_qc_show_imputed == true || params.fp_analyst_qc_show_imputed == 'true') ? 'true' : 'false'
    def qc_both = (params.fp_analyst_qc_include_both == true || params.fp_analyst_qc_include_both == 'true') ? 'true' : 'false'
    def volcano_display_names = (params.fp_analyst_volcano_display_names == true || params.fp_analyst_volcano_display_names == 'true') ? 'true' : 'false'
    def volcano_show_gene = (params.fp_analyst_volcano_show_gene == true || params.fp_analyst_volcano_show_gene == 'true') ? 'true' : 'false'
    def volcano_highlight_feature = (params.fp_analyst_volcano_highlight_feature != null && params.fp_analyst_volcano_highlight_feature != '') ? params.fp_analyst_volcano_highlight_feature : ''
    """
    # Create output directory
    mkdir -p output/

    # Run FragPipe-Analyst R script (single run with both protein and gene feature lists)
    Rscript "${script_path}" \\
        --experiment_annotation "${experiment_annotation}" \\
        --quantification_file "${quantification_file}" \\
        --mode "${mode}" \\
        --level "${data_type}" \\
        --feature_list_protein "${feature_list_protein}" \\
        --feature_list_gene "${feature_list_gene}" \\
        --top_n_protein "${top_n_protein}" \\
        --top_n_gene "${top_n_gene}" \\
        --pathway_database "${pathway_database}" \\
        --pathway_direction "${pathway_direction}" \\
        --go_database "${go_database}" \\
        --go_direction "${go_direction}" \\
        --lfq_type "${lfq_type}" \\
        --de_alpha "${de_alpha}" \\
        --de_lfc "${de_lfc}" \\
        --de_fdr "${de_fdr}" \\
        --min_global_appearance "${min_global}" \\
        --min_appearance_one_condition "${min_cond}" \\
        --qc_show_imputed "${qc_imputed}" \\
        --qc_include_both "${qc_both}" \\
        --volcano_display_names "${volcano_display_names}" \\
        --volcano_show_gene "${volcano_show_gene}" \\
        --volcano_highlight_feature "${volcano_highlight_feature}" \\
        --output_dir "output/"
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    FragPipeAnalystR: \$(Rscript -e 'tryCatch({cat(as.character(packageVersion(\"FragPipeAnalystR\")))}, error=function(e){cat(\"unknown\")})' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}
