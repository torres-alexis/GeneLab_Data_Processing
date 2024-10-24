process FETCH_ISA {

    tag "OSD-${osd}"

    publishDir params.glds ? "${params.outdir}/GLDS-${params.glds}/Metadata" : "${params.outdir}/Metadata",
        mode: params.publish_dir_mode

    input:
    val(osd)
    val(glds)

    output:
    path "*.zip", emit: isa_archive

    script:
    """
    python $projectDir/bin/fetch_isa.py --osd ${osd} --outdir .
    """
}
