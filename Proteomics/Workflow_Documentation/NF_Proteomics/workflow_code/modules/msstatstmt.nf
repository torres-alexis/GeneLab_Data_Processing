process MSSTATSTMT {
    publishDir path: { "${output_dir}/MSstatsTMT/" },
        mode: params.publish_dir_mode,
        pattern: "msstatstmt_comparison*.csv"
    publishDir path: { "${output_dir}/MSstatsTMT/" },
        mode: params.publish_dir_mode,
        pattern: "msstatstmt_contrasts*.csv"
    publishDir path: { "${output_dir}/MSstatsTMT/" },
        mode: params.publish_dir_mode,
        pattern: "dropped-conditions-msstatstmt*.txt"
    publishDir path: { "${output_dir}/MSstatsTMT/" },
        mode: params.publish_dir_mode,
        pattern: "dropped-runs-msstatstmt*.txt"

    input:
    val(output_dir)
    path(msstats_tmt_annotation)
    path(msstats_csv)

    output:
    path("versions.yml"), emit: versions
    path("msstatstmt_comparison*.csv"), emit: comparison, optional: true
    path("msstatstmt_contrasts*.csv"), emit: contrasts, optional: true
    path("dropped-conditions-msstatstmt*.txt"), emit: conditions_notice, optional: true
    path("dropped-runs-msstatstmt*.txt"), emit: runs_notice, optional: true

    script:
    """
    msstatstmt_analysis.R . ${msstats_tmt_annotation} ${msstats_csv} ${params.assay_suffix}

    echo '"${task.process}":' > versions.yml
    echo "    MSstatsTMT: \$(Rscript -e 'cat(as.character(packageVersion(\"MSstatsTMT\")))' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}
