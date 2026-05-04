process PACKAGE_PROCESSING_INFO { 
    publishDir "${ch_outdir}/GeneLab",
        mode: params.publish_dir_mode,
        pattern: "*.zip"

    input:
        path(processing_scripts)
        val(ch_outdir)

    output:
        path("processing_info${params.assay_suffix}.zip"), emit: zip

    script:
    """
    for f in ${processing_scripts}/nextflow*.txt; do
        echo "Purging file paths from \$f"
        clean_paths.sh "\$f"

        # Add assay suffix to the end of nextflow log files before .txt extension, if not already present
        if [[ "\$f" != *"${params.assay_suffix}.txt" ]]; then
            mv "\$f" "\${f%.txt}${params.assay_suffix}.txt"
        fi
    done
        
    # Zip 
    zip -r processing_info${params.assay_suffix}.zip processing_scripts
    """
}