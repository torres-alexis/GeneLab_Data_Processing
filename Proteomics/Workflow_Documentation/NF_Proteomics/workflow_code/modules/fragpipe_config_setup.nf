process FRAGPIPE_CONFIG_SETUP {
    tag "${workflow_config.getName()}"

    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "output/*.workflow",
        saveAs: { it.replace("output/", "") }

    input:
    val(ch_outdir)
    path(workflow_config)
    path(proteome)

    output:
    path("output/*.workflow"), emit: fragpipe_config
    stdout emit: fragpipe_json

    script:
    def assay_suffix_flag = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    """
    fragpipe_config_setup.py \\
        --input ${workflow_config} \\
        --output output/ \\
        --proteome ${proteome} \\
        ${assay_suffix_flag} \\
        --stdout
    """
}