/*
 * Different Gene Expression Analysis Processes
 */
// ERCC counts are removed 

process DESEQ2_DGE {

    input:
        path(runsheet_path)
        path(gene_counts)
        val(meta)

    output: 
        tuple path("Normalized_Counts_GLbulkRNAseq.csv"),
                path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts_GLbulkRNAseq.csv" : 
                "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"),         emit: norm_counts
        tuple path("contrasts_GLbulkRNAseq.csv"),
                path("SampleTable_GLbulkRNAseq.csv"),
                path("differential_expression_GLbulkRNAseq.csv"),
                path("visualization_output_table_GLbulkRNAseq.csv"),
                path("visualization_PCA_table_GLbulkRNAseq.csv"),       emit: dge
                path("summary.txt"),                                    emit: summary
                path("versions.txt"),                                   emit: versions_txt

    script:
        """
        Rscript --vanilla deseq2_dge.R \\
            --runsheet_path ${ runsheet_path } \\
            --gene_id_type ${ meta.gene_id_type } \\
            ${ params.mode == 'microbes' ? '--microbes' : ''} \\
            ${ params.use_dummy_gene_counts ? '--DEBUG_MODE_ADD_DUMMY_COUNTS' : ''} \\
            --input_gene_results_dir ${ gene_counts } \\
            --normalization 'default' 
        """
}