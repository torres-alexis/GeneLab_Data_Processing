process RUNSHEET_TO_MANIFEST {
    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "*.tsv"

    input: 
    val(ch_outdir)
    path(runsheet)

    output:
    path("*.tsv"), emit: manifest

    script:
    def assay_suffix_flag = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
    def mode = (params.fragpipe_workflow && params.fragpipe_workflow.startsWith('TMT')) ? "TMT" : "LFQ"
    def mode_flag = "--mode ${mode}"
    def output_filename = "manifest${params.assay_suffix ?: ''}.tsv"
    """
    runsheet_to_manifest.py --runsheet ${runsheet} --output ${output_filename} ${assay_suffix_flag} ${mode_flag}
    """
}