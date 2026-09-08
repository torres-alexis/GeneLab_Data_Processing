process ISA_TO_RUNSHEET {
    tag "${osd_accession}_${glds_accession}"

    input: 
    val(ch_outdir)
    val(osd_accession)
    val(glds_accession)
    path(isa_archive)
    path(dp_tools_plugin)

    output:
    path("*.csv"), emit: runsheet
    path("isa_archive/${isa_archive}"), emit: isa_copy
    //path("versions.yml"), emit: versions

    script:
    """
    dpt-isa-to-runsheet --accession ${glds_accession} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}

    # Copy the ISA archive to the output directory
    mkdir -p isa_archive
    cp ${isa_archive} isa_archive/

    #echo '"${task.process}":' > versions.yml
    #echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}