process VALIDATE_PROCESSING {
    publishDir "${ch_outdir}/GeneLab",
        mode: params.publish_dir_mode,
        pattern: "*.log"

    input:
        path(ch_outdir)
        path(raw_md5sum)
        path(processed_md5sum)

    output:
        path("validate_processed_proteomics${params.assay_suffix}.log"), emit: validation_log

    script:
        """
        GL-validate-processed-proteomics-data.py \\
            --outdir ${ch_outdir} \\
            --assay_suffix ${params.assay_suffix} \\
            --output validate_processed_proteomics${params.assay_suffix}.log
        """
}
