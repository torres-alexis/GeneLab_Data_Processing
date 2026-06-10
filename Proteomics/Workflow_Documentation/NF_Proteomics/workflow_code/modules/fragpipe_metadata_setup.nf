process FRAGPIPE_METADATA_SETUP {
    publishDir path: { "${ch_outdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "*.tsv"
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
    path("msstats_tmt_annotation*.tsv"), emit: msstats_tmt_annotation, optional: true
    path("versions.yml"), emit: versions
    path("sheets/*"), emit: sheets

    script:
    def is_tmt = params.fragpipe_workflow?.contains('TMT')
    def assay_suffix_flag = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    def sheet_flag = is_tmt ? "--data_sheet ${sheets[0]} --sample_sheet ${sheets[1]}" : "--runsheet ${sheets[0]}"
    def msstats_anno_out = params.assay_suffix ?
        "msstats_tmt_annotation${params.assay_suffix}.tsv" :
        "msstats_tmt_annotation.tsv"
    """
    runsheet_to_fp_metadata.py ${sheet_flag} ${assay_suffix_flag}

    if [[ "${params.fragpipe_workflow}" == TMT* ]]; then
      build_msstats_tmt_annotation.py \\
        --data_sheet ${sheets[0]} \\
        --sample_sheet ${sheets[1]} \\
        --output ${msstats_anno_out}
    fi

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
