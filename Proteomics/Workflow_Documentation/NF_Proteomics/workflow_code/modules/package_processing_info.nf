process PACKAGE_PROCESSING_INFO { 
    publishDir path: { "${ch_outdir}/GeneLab" },
        mode: params.publish_dir_mode,
        pattern: "*.zip"

    input:
        path(processing_info)
        val(ch_outdir)

    output:
        path("processing_info${params.assay_suffix}.zip"), emit: zip

    script:
    """
    for f in ${processing_info}/nextflow*.txt; do
        # Add assay suffix to the end of nextflow log files before .txt extension, if not already present
        if [[ "\$f" != *"${params.assay_suffix}.txt" ]]; then
            mv "\$f" "\${f%.txt}${params.assay_suffix}.txt"
        fi
    done
        
    # Zip 
    zip -r processing_info${params.assay_suffix}.zip processing_info
    """
}