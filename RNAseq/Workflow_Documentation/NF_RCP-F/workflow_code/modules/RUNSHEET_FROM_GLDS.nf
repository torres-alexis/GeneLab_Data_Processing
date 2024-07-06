process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.{zip,csv}",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path(dp_tools_plugin)

  output:
    path("${ osdAccession }_*_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    """

    dpt-get-isa-archive --accession ${ osdAccession }
    ls ${dp_tools_plugin}

    dpt-isa-to-runsheet --accession ${ osdAccession } \
      --plugin-dir ${dp_tools_plugin} \
      --isa-archive *.zip ${ injects }
    """
}