process MOVE_RUNSHEET {
    tag "${ gldsAccession }"
    publishDir "${params.outputDir}/${params.gldsAccession}/Metadata",
        pattern: "*.csv",
        mode: params.publish_dir_mode

    input:
        path(runsheetPath)

    output:
        path(runsheetPath), emit: runsheet

    script:
        """
        # Move the provided runsheet to the specified publish directory
        touch ${ runsheetPath }
        """
}