process VV_STEP {
    tag "${name}"

    publishDir path: { "${output_dir}/VV_Logs" },
        mode: params.publish_dir_mode,
        pattern: "VV_log_*.log"

    input:
        val output_dir
        val name
        path files

    output:
        path("VV_log_*.log"), emit: log

    script:
    def suffix = (params.assay_suffix != null && params.assay_suffix != '') ? params.assay_suffix.toString() : ''
    """
    vv_step.py --name ${name} --suffix "${suffix}" ${files}
    """
}
