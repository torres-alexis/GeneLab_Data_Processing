process MSSTATS_TMT {
    publishDir path: { "${output_dir}/MSstats/" },
        mode: params.publish_dir_mode,
        pattern: "msstats_comparison*.csv"
    publishDir path: { "${output_dir}/MSstats/" },
        mode: params.publish_dir_mode,
        pattern: "msstats_contrasts*.csv"

    input:
    val(output_dir)
    path(msstats_tmt_annotation)
    path(msstats_tmt_csv)

    output:
    path("versions.yml"), emit: versions
    path("msstats_comparison*.csv"), emit: comparison, optional: true
    path("msstats_contrasts*.csv"), emit: contrasts, optional: true

    script:
    """
    mkdir -p msstats_inputs
    for f in ${msstats_tmt_csv}; do
      base=\$(basename "\$(dirname "\$f")")
      cp "\$f" "msstats_inputs/\${base}_msstats.csv"
    done

    msstats_tmt_analysis.R . ${msstats_tmt_annotation} msstats_inputs ${params.assay_suffix}

    echo '"${task.process}":' > versions.yml
    echo "    msstatstmt: \$(Rscript -e 'cat(as.character(packageVersion(\"MSstatsTMT\")))' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}
