process CLEAN_FRAGPIPE_TABLES {
    tag "fragpipe_tables"

    input:
    val(output_dir)
    path(fragpipe_tables)
    path(experiment_annotation)

    output:
    path("output/*.tsv"), emit: tables, optional: true

    script:
    def fragpipe_table_args = fragpipe_tables instanceof List ? fragpipe_tables.join(' ') : fragpipe_tables
    def analysis_type = params.fragpipe_workflow?.startsWith('TMT') ? 'TMT' : 'LFQ'
    """
    python ${projectDir}/bin/clean_fragpipe_tables.py \\
        --experiment_annotation ${experiment_annotation} \\
        --output_dir output \\
        --type ${analysis_type} \\
        ${fragpipe_table_args}
    """
}
