process FETCH_ISA {

    tag "${osd_accession}"

    input:
    val(output_dir)
    val(osd_accession)
    val(glds_accession)
    output:
    path "*.zip", emit: isa_archive

    script:
    """
    fetch_isa.py --osd ${osd_accession} --outdir .
    """
}