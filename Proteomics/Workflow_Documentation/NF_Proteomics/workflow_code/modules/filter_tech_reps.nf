process FILTER_TECH_REPS {
    tag "${sheet_file.name}"

    input:
        val ch_outdir
        path sheet_file

    output:
        path("original/${sheet_file.name}"), emit: original
        path('filtered_sheet.csv'), emit: filtered
        path('publish/filtered_sheet.csv'), emit: filtered_publish, optional: true
        path('tech_reps_dropped.tsv'), emit: drop_log, optional: true

    script:
    def mode = params.fragpipe_workflow?.startsWith('TMT') ? 'tmt' : 'lfq'
    """
    mkdir -p original
    cp "${sheet_file}" "original/${sheet_file.name}"

    filter_tech_reps.py \\
        --mode ${mode} \\
        --input ${sheet_file} \\
        --output filtered_sheet.csv \\
        --drop-log tech_reps_dropped.tsv \\
        --publish-if-changed publish/filtered_sheet.csv
    """
}
