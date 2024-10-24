process MULTIQC {
    tag("Dataset-wide")
    
    input:
    path(sample_names)
    path(fastqc_zips, stageAs: 'fastqc/*')
    path(multiqc_config)
    

    output:
    // path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip"), emit: zipped_report
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report"), emit: unzipped_report
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq.html"), emit: html
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data"), emit: data
    
    
    path("versions.yml"), emit: versions

    script:
    def config_arg = multiqc_config.name != "NO_FILE" ? "--config ${multiqc_config}" : ""
    """
    multiqc \\
        --force \\
        --interactive \\
        -o ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report \\
        -n ${ params.MQCLabel }_multiqc_GLbulkRNAseq \\
        $config_arg \\
        .
    # zip -r '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip' '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report'
    # Use awk to clean paths in relevant multiqc output files, leaving them starting from the nextflow '/work/' directory
    awk '{gsub(/\\/.*\\/work\\//, "/work/"); print}' ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc_data.json > temp_data.json && mv temp_data.json ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc_data.json
    awk '{gsub(/\\/.*\\/work\\//, "/work/"); print}' ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc_sources.txt > temp_sources.txt && mv temp_sources.txt ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc_sources.txt
    awk '{gsub(/\\/.*\\/work\\//, "/work/"); print}' ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc.log > temp_log.txt && mv temp_log.txt ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data/multiqc.log

    echo "${task.process}:" > versions.yml
    echo "    multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")" >> versions.yml
    """
}
