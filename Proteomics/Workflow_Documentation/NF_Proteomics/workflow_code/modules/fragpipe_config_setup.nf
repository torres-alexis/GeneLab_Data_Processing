process FRAGPIPE_CONFIG_SETUP {
    tag "${workflow_config.getName()}"

    input:
    val(ch_outdir)
    path(workflow_config)
    path(proteome)

    output:
    path("output/*.workflow"), emit: fragpipe_config
    stdout emit: fragpipe_json

    script:
    def assay_suffix_flag = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    def tmt_flag = params.fragpipe_workflow?.startsWith('TMT') ? '--tmt' : ''
    def tmt_extraction_tool_flag = params.fragpipe_workflow?.startsWith('TMT') ?
        "--tmt_extraction_tool ${params.tmt_extraction_tool}" : ''
    """
    fragpipe_config_setup.py \\
        --input ${workflow_config} \\
        --output output/ \\
        --proteome ${proteome} \\
        ${assay_suffix_flag} \\
        ${tmt_flag} \\
        ${tmt_extraction_tool_flag} \\
        --stdout
    """
}