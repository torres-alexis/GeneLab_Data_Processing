process ISA_TO_RUNSHEET {
    tag "${osd_accession}_${glds_accession}"

    publishDir "${glds_accession}/Metadata",
        mode: params.publish_dir_mode

    input: 
    val(osd_accession)
    val(glds_accession)
    path isa_archive
    path dp_tools_plugin

    output:
    path "*.csv", emit: runsheet

    script:
    """
    dpt-isa-to-runsheet --accession ${osd_accession} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}
    """
}