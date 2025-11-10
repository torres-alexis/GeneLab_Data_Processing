process RAWBEANS_QC {
    tag "${meta.id}"
    
    publishDir "${output_dir}/00-RawData/${meta.id}/",
        mode: params.publish_dir_mode,
        pattern: "${meta.id}${params.assay_suffix}/qc-report.html",
        saveAs: { filename -> "${meta.id}${params.assay_suffix}_qc-report.html" }
    publishDir "${output_dir}/00-RawData/${meta.id}/",
        mode: params.publish_dir_mode,
        pattern: "${meta.id}${params.assay_suffix}/resources/**",
        saveAs: { filename -> "resources/**" }

    input:
    val(output_dir)
    tuple val(meta), path(mzml_file)

    output:
    tuple val(meta), path("${meta.id}${params.assay_suffix}/qc-report.html"), emit: qc_report
    path("${meta.id}${params.assay_suffix}/resources/**"), emit: resources
    path("versions.yml"), emit: versions

    script:
    """
    # Run RawBeans QC
    create-qc-report.py \\
        --input ${mzml_file} \\
        --output-dir . \\
        --batch \\
        --cores ${task.cpus}

    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    protqc: \$(pip show protqc | grep Version | cut -d' ' -f2)" >> versions.yml
    echo "    python: \$(python --version 2>&1 | sed 's/Python //g')" >> versions.yml
    """
}

process RAWBEANS_QC_ALL {
    
    publishDir "${output_dir}/00-RawData/All-rawbeans/",
        mode: params.publish_dir_mode,
        pattern: "qc-report.html",
        saveAs: { filename -> "All-rawbeans_qc-report.html" }

    publishDir "${output_dir}/00-RawData/All-rawbeans/",
        mode: params.publish_dir_mode,
        pattern: "resources/**"
    
    input:
    val(output_dir)
    path(mzml_files)

    output:
    path("qc-report.html"), emit: qc_report
    path("resources/**"), emit: resources
    path("versions.yml"), emit: versions

    script:
    """
    # Run RawBeans QC on all samples
    create-qc-report.py \\
        --input ${mzml_files} \\
        --output-dir . \\
        --cores ${task.cpus}

    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    protqc: \$(pip show protqc | grep Version | cut -d' ' -f2)" >> versions.yml
    echo "    python: \$(python --version 2>&1 | sed 's/Python //g')" >> versions.yml
    """
}