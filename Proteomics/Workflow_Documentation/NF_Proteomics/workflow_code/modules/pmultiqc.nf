process PMULTIQC {
    publishDir path: { "${publishdir}" },
        pattern:  "*.{html,zip}" ,
        mode: params.publish_dir_mode
    
    input:
    val(publishdir)
    path(fragpipe_output_dir)

    output:
    path("multiqc${ params.assay_suffix }.html"), emit: html
    path("multiqc${ params.assay_suffix }_data.zip"), emit: zipped_data

    path("versions.yml"), emit: versions

    script:
    """
    multiqc \
        --fragpipe-plugin \
        -c ${params.multiqc_config} \
        -o . \
        -n multiqc${ params.assay_suffix } \
        -z \
        "${fragpipe_output_dir}"
    
    # Create versions.yml
    echo '"${task.process}":' > versions.yml
    echo '    multiqc: '\$(/usr/local/bin/multiqc --version 2>&1 | sed 's/multiqc, version //' || echo 'unknown') >> versions.yml
    echo '    pmultiqc: '\$(/usr/local/bin/multiqc --pmultiqc-version 2>&1 | sed 's/pmultiqc, version //' || echo 'unknown') >> versions.yml
    """
}
