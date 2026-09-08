process VALIDATE_PROCESSING {
    publishDir path: { "${ch_outdir}/GeneLab" },
        mode: params.publish_dir_mode,
        pattern: "*.log"

    input:
        val(ch_outdir)
        path(raw_md5sum)
        path(processed_md5sum)

    output:
        path("validate_processed_proteomics${params.assay_suffix}.log"), emit: validation_log

    script:
        def decoy_prefix = (params.philosopher_decoy_prefix != null && params.philosopher_decoy_prefix != '') ? params.philosopher_decoy_prefix : 'rev_'
        """
        GL-validate-processed-proteomics-data.py \\
            --outdir ${ch_outdir} \\
            --assay_suffix ${params.assay_suffix} \\
            --decoy-prefix "${decoy_prefix}" \\
            --output validate_processed_proteomics${params.assay_suffix}.log
        """
}
