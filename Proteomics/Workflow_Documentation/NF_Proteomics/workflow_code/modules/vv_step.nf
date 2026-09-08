process VV_STEP {
    tag "${name}"

    input:
        val output_dir
        tuple val(name), path(files)

    output:
        path("VV_log_*.log"), emit: log

    script:
    def suffix = (params.assay_suffix != null && params.assay_suffix != '') ? params.assay_suffix.toString() : ''
    """
    vv_step.py --name ${name} --suffix "${suffix}" ${files}
    """
}
