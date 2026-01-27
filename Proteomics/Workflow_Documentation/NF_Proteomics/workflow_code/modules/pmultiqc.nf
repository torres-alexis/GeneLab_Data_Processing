process PMULTIQC {
    // tag("Dataset-wide")
    // publishDir "${ publishdir }",
    //     pattern:  "*.{html,zip}" ,
    //     mode: params.publish_dir_mode
    
    input:
    val(publishdir)
    path(fragpipe_output_dir) // FragPipe output directory containing psm.tsv files
    

    output:
    path("multiqc${ params.assay_suffix }.html"), emit: html
    path("multiqc${ params.assay_suffix }_data"), emit: data
    path("multiqc${ params.assay_suffix }_data.zip"), emit: zipped_data

    path("versions.yml"), emit: versions

    script:
    """
    # Unset any config-related environment variables that might interfere
    unset MULTIQC_CONFIG || true
    
    # Call multiqc directly (entrypoint is overridden to /bin/sh)
    /usr/local/bin/multiqc \
        --fragpipe-plugin \
        -o . \
        -n multiqc${ params.assay_suffix } \
        "${fragpipe_output_dir}"
    
    # Clean paths and create zip
    clean_multiqc_paths.py multiqc${ params.assay_suffix }_data .
    
    # Create versions.yml
    echo '"${task.process}":' > versions.yml
    echo '    multiqc: '\$(/usr/local/bin/multiqc --version 2>&1 | sed 's/multiqc, version //' || echo 'unknown') >> versions.yml
    echo '    pmultiqc: '\$(/usr/local/bin/multiqc --pmultiqc-version 2>&1 | sed 's/pmultiqc, version //' || echo 'unknown') >> versions.yml
    """
}
