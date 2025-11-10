process MSSTATS {
    publishDir "${output_dir}/05-MSstats/",
        mode: params.publish_dir_mode,
        pattern: "msstats_input.csv",
        saveAs: { filename -> "msstats.csv" }
    publishDir "${output_dir}/05-MSstats/",
        mode: params.publish_dir_mode,
        pattern: "msstats_comparison_*.csv"
    publishDir "${output_dir}/05-MSstats/",
        mode: params.publish_dir_mode,
        pattern: "msstats_comparison_all.csv"
    publishDir "${output_dir}/05-MSstats/",
        mode: params.publish_dir_mode,
        pattern: "msstats_contrasts.csv"

    input:
    val(output_dir)
    path(msstats_csv)
    path(runsheet)

    output:
    path("versions.yml"), emit: versions
    path("msstats_input.csv"), emit: msstats_processed
    path("msstats_comparison_*.csv"), emit: comparison, optional: true
    path("msstats_comparison_all.csv"), emit: comparison_all, optional: true
    path("msstats_contrasts.csv"), emit: contrasts, optional: true

    script:
    """
    msstats_analysis.R . ${params.assay_suffix ?: ''} ${runsheet} ${msstats_csv}
    
    # Version info (back in work directory)
    echo '"${task.process}":' > versions.yml
    echo "    msstats: \$(Rscript -e 'cat(as.character(packageVersion(\"MSstats\")))' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}

