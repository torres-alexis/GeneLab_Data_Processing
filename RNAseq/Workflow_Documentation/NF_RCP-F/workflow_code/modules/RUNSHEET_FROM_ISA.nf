process RUNSHEET_FROM_ISA {
  // Generates Runsheet using a path to an ISA archive
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.{zip,csv}",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path(isaArchivePath)
    path(dp_tools_plugin)

  output:
    path("${ osdAccession }_*_v?_runsheet.csv"), emit: runsheet
    path(isaArchivePath), emit: isaArchive

  script:
    """
    dpt-isa-to-runsheet --accession ${ osdAccession } \
      --plugin-dir ${dp_tools_plugin} \
      --isa-archive ${ isaArchivePath }
    """
}