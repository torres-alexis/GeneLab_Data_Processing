process ISA_TO_RUNSHEET {
    tag "OSD-${params.osd}_GLDS-${glds}"

    publishDir params.glds ? "${params.outdir}/GLDS-${params.glds}/Metadata" : "${params.outdir}/Metadata",
        mode: params.publish_dir_mode

    input: 
    val(glds)
    path isa_archive
    path dp_tools_plugin

    output:
    path "*.csv", emit: runsheet

    script:
    """
    dpt-isa-to-runsheet --accession GLDS-${glds} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}
    """
}