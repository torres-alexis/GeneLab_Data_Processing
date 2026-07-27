process FRAGPIPE_METADATA_SETUP {
    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "*.tsv"
    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "MSstatsTMT_annotation*.csv"
    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "sheets/*",
        saveAs: { filename ->
            if (filename.startsWith("sheets/")) return filename.replace("sheets/", "")
            else return filename
        }

    input:
    val(ch_outdir)
    path(sheets)

    output:
    path("manifest*.tsv"), emit: manifest
    path("experiment_annotation*.tsv"), emit: experiment_annotation, optional: true
    path("MSstatsTMT_annotation*.csv"), emit: msstats_tmt_annotation, optional: true
    path("versions.yml"), emit: versions
    path("sheets/*"), emit: sheets

    script:
    def assay_suffix_flag = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    def workflow_flag = params.fragpipe_workflow ? "--fragpipe_workflow ${params.fragpipe_workflow}" : ""
    def sheet_flag = params.fragpipe_workflow?.contains('TMT') ? "--data_sheet ${sheets[0]} --sample_sheet ${sheets[1]}" : "--runsheet ${sheets[0]}"
    """
    runsheet_to_metadata.py ${sheet_flag} ${assay_suffix_flag} ${workflow_flag}

    # Create output dir and copy input sheet(s) there for publishing
    mkdir -p sheets
    cp "${sheets[0]}" sheets/
    if [[ "${params.fragpipe_workflow}" == TMT* ]]; then
      cp "${sheets[1]}" sheets/
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools 2>/dev/null | grep '^Version:' | sed 's/Version: //')" >> versions.yml
    """
}
