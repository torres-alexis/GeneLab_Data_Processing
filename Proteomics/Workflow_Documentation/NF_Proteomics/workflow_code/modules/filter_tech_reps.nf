process FILTER_TECH_REPS {
    tag "${sheet_file.name}"

    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "original/*",
        saveAs: { fname -> fname.replaceFirst(/^original\//, '') }

    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "publish/*",
        saveAs: { fname -> fname.replaceFirst(/^publish\//, '') }

    input:
        val ch_outdir
        path sheet_file

    output:
        path("original/${sheet_file.name}"), emit: original
        path('filtered_sheet.csv'), emit: filtered
        path('publish/filtered_sheet.csv'), emit: filtered_publish, optional: true

    script:
    def mode = params.fragpipe_workflow?.startsWith('TMT') ? 'tmt' : 'lfq'
    """
    mkdir -p original
    cp "${sheet_file}" "original/${sheet_file.name}"

    filter_tech_reps.py \\
        --mode ${mode} \\
        --input ${sheet_file} \\
        --output filtered_sheet.csv \\
        --publish-if-changed publish/filtered_sheet.csv
    """
}
