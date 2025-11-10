process OPENMS_FILEINFO {
    tag "${meta.id}"
    
    publishDir "${output_dir}/01-QC/${meta.id}",
        mode: params.publish_dir_mode,
        pattern: "*.txt"

    input:
    val(output_dir)
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("*_fileinfo.txt"), emit: fileinfo_report
    path("versions.yml"), emit: versions

    script:
    """
    # Run OpenMS FileInfo for mzML validation and metadata extraction
    FileInfo -in ${mzml_file} -out ${meta.id}_fileinfo.txt
    
    # Check if FileInfo succeeded (validates mzML structure)
    if [ \$? -eq 0 ]; then
        echo "mzML validation: PASSED" >> ${meta.id}_fileinfo.txt
        echo "File: ${mzml_file}" >> ${meta.id}_fileinfo.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_fileinfo.txt
    else
        echo "mzML validation: FAILED" >> ${meta.id}_fileinfo.txt
        echo "ERROR: FileInfo failed for ${meta.id}" >&2
        exit 1
    fi
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    openms: \$(FileInfo 2>&1 | grep -o 'Version [0-9.]*' | head -1 || echo 'unknown')" >> versions.yml
    """
}
