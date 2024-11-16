/*
 * Different Gene Expression Analysis Processes
 */
// ERCC counts are removed before normalization

process DESEQ2_DGE {
    tag "Dataset-wide"

    input:
        path(runsheet_path)
        path(gene_counts)
        val(meta)

    output:
        tuple path("Normalized_Counts${params.output_suffix}.csv"),
              path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts${params.output_suffix}.csv" : 
                   "RSEM_Unnormalized_Counts${params.output_suffix}.csv"),                emit: norm_counts
        tuple path("contrasts${params.output_suffix}.csv"),
              path("SampleTable${params.output_suffix}.csv"),
              path("differential_expression_no_annotations${params.output_suffix}.csv"),  emit: dge
        path("VST_Normalized_Counts${params.output_suffix}.csv"),                         emit: vst_norm_counts
        path("summary.txt"),                                                              emit: summary
        path("DESeq2_DGE.html"),                                                          emit: dge_html
        path("versions.txt"),                                                             emit: versions_txt

    script:
        def dge_rmd_file = "${projectDir}/bin/deseq2_dge.Rmd"
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE' : 'FALSE'
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'
        def output_filename_suffix = params.output_suffix ?: ""

        """
        Rscript -e "rmarkdown::render('${dge_rmd_file}', 
        output_file = 'DESeq2_DGE.html',
        output_dir = '\${PWD}',
            params = list(
                cpus = ${task.cpus},
                runsheet_path = '${runsheet_path}',
                input_gene_results_dir = '\${PWD}',
                gene_id_type = '${meta.gene_id_type}',
                microbes = ${microbes},
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts}
            ))"
        """
}