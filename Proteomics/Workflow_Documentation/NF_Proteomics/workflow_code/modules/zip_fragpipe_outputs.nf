process ZIP_FRAGPIPE_OUTPUTS {
    tag "fragpipe_zip"

    publishDir path: { "${output_dir}/FragPipe/" },
        mode: params.publish_dir_mode,
        pattern: "*.zip"

    input:
    val(output_dir)
    path(fragpipe_output_dir)

    output:
    path("*.zip"), emit: zip, optional: true

    script:
    def workflow_type = params.fragpipe_workflow?.startsWith('TMT') ? 'TMT' : 'LFQ'
    def zip_name = "fragpipe${params.assay_suffix ?: ''}.zip"
    """
    python ${projectDir}/bin/zip_fragpipe_outputs.py \\
        --input ${fragpipe_output_dir} \\
        --output ${zip_name} \\
        --workflow-type ${workflow_type}
    """
}
