process VV_STEP {
    tag "${name}"

    input:
        val output_dir
        tuple val(name), path(files)

    output:
        path("VV_log_*.csv"), emit: log

    script:
    def suffix = (params.assay_suffix != null && params.assay_suffix != '') ? params.assay_suffix.toString() : ''
    """
    vv_step.py --name ${name} --suffix "${suffix}" ${files}
    """
}

process VV_CONCAT_FILTER {
    input:
        val output_dir
        path("VV_in.csv")

    output:
        tuple path("VV_log_final${params.assay_suffix}.csv"), path("VV_log_final_only_issues${params.assay_suffix}.csv"), emit: logs

    script:
    """
    concat_logs.py --assay_suffix ${params.assay_suffix}
    filter_to_only_issues.py --assay_suffix ${params.assay_suffix}
    """
}
