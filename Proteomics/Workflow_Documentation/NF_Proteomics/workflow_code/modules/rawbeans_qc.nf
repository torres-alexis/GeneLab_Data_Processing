process RAWBEANS_QC {
    tag "${meta.id}"
    
    publishDir "${output_dir}/RawBeans_QC/${meta.id}/",
        mode: params.publish_dir_mode,
        pattern: "${meta.id}${params.assay_suffix}_qc-report.zip"

    input:
    val(output_dir)
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("${meta.id}${params.assay_suffix}_qc-report.zip"), emit: qc_report
    path("versions.yml"), emit: versions

    script:
    """
    # Run RawBeans QC
    create-qc-report.py \\
        --input ${mzml_file} \\
        --output-dir . \\
        --batch \\
        --cores ${task.cpus}

    # Create zip file with HTML and resources folder (cd to get files at root level)
    cd ${meta.id}${params.assay_suffix}
    zip -r ../${meta.id}${params.assay_suffix}_qc-report.zip qc-report.html resources/
    cd ..

    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    protqc: \$(pip show protqc | grep Version | cut -d' ' -f2)" >> versions.yml
    echo "    python: \$(python --version 2>&1 | sed 's/Python //g')" >> versions.yml
    """
}

process RAWBEANS_QC_ALL {
    
    publishDir "${output_dir}/RawBeans_QC/",
        mode: params.publish_dir_mode,
        pattern: "All${params.assay_suffix}_qc-report.zip"
    
    input:
    val(output_dir)
    path(mzml_files)

    output:
    path("All${params.assay_suffix}_qc-report.zip"), emit: qc_report
    path("versions.yml"), emit: versions

    script:
    """
    # Run RawBeans QC on all samples
    create-qc-report.py \\
        --input ${mzml_files} \\
        --output-dir . \\
        --cores ${task.cpus}

    # Create zip file with HTML and resources folder
    zip -r All${params.assay_suffix}_qc-report.zip qc-report.html resources/

    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    protqc: \$(pip show protqc | grep Version | cut -d' ' -f2)" >> versions.yml
    echo "    python: \$(python --version 2>&1 | sed 's/Python //g')" >> versions.yml
    """
}