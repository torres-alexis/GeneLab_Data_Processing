process RUNSHEET_TO_EXPERIMENT_ANNOTATION {
    tag "LFQ"
    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "experiment_annotation.tsv"

    input:
    val(ch_outdir)
    path(experiment_annotation)
    path(runsheet)

    output:
    path("experiment_annotation.tsv"), emit: experiment_annotation

    script:
    """
    ${projectDir}/bin/runsheet_to_experiment_annotation.py \\
        --experiment_annotation ${experiment_annotation} \\
        --runsheet ${runsheet} \\
        --output experiment_annotation.tsv

    """
}
