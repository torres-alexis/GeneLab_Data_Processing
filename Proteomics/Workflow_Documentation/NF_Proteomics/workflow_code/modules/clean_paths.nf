process CLEAN_PATHS {
    input:
        val(processed_dir)

    output:
        val(processed_dir), emit: processed_dir

    script:
    """
    python ${projectDir}/bin/clean_paths.py \\
        --processed-dir "${processed_dir}" \\
        --assay-suffix "${params.assay_suffix}"
    """
}
