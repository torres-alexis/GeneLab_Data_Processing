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
        tuple path("Normalized_Counts_GLbulkRNAseq.csv"),
             path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts_GLbulkRNAseq.csv" : 
                  "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),    emit: norm_counts
        tuple path("contrasts_GLbulkRNAseq.csv"),
             path("SampleTable_GLbulkRNAseq.csv"),
             path("differential_expression_no_annotations_GLbulkRNAseq.csv"), emit: dge
        path("summary.txt"),                                       emit: summary
        path("versions.txt"),                                      emit: versions_txt
        path("DESeq2_DGE.html"),                                   emit: dge_html

    script:
        def dge_rmd_file = "${projectDir}/bin/deseq2_dge.Rmd"
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE' : 'FALSE'
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'

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
                normalization = 'default',
                work_dir = '\${PWD}',
                SUMMARY_FILE_PATH = 'summary.txt',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts}
            ))"
        """
}