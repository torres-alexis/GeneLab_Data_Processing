process MSSTATS {
    // MSstats inputs (e.g. msstats.csv, msstats_ptm.csv) are published under FragPipe/
    // This process only publishes MSstats results here.
    publishDir path: { "${output_dir}/MSstats/" },
        mode: params.publish_dir_mode,
        pattern: "msstats_comparison*.csv"
    publishDir path: { "${output_dir}/MSstats/" },
        mode: params.publish_dir_mode,
        pattern: "msstats_contrasts*.csv"

    input:
    val(output_dir)
    path(experiment_annotation)
    path(msstats_csv)

    output:
    path("versions.yml"), emit: versions
    path("msstats_comparison*.csv"), emit: comparison, optional: true
    path("msstats_contrasts*.csv"), emit: contrasts, optional: true

    script:
    def drop_dc = (params.drop_decoys_contams != false && params.drop_decoys_contams != 'false') ? 'true' : 'false'
    def decoy_prefix = (params.philosopher_decoy_prefix != null && params.philosopher_decoy_prefix != '') ? params.philosopher_decoy_prefix : 'rev_'
    """
    export DROP_DECOYS_CONTAMS=${drop_dc}
    export PHILOSOPHER_DECOY_PREFIX=${decoy_prefix}
    msstats_analysis.R . ${experiment_annotation} ${msstats_csv} ${params.assay_suffix}
    
    # Version info (back in work directory)
    echo '"${task.process}":' > versions.yml
    echo "    msstats: \$(Rscript -e 'cat(as.character(packageVersion(\"MSstats\")))' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}

