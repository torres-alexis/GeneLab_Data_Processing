process CONCAT_ERCC {
    // Concanates ERCC fasta and gtf to reference fasta and gtf
    storeDir "${reference_store_path}/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"

    input:
        val(reference_store_path)
        val(organism_sci)
        val(reference_source)
        val(reference_version)
        path(genome_fasta)
        path(genome_gtf)
        path(ercc_fasta)
        path(ercc_gtf)
        val(has_ercc)


    output:
        path("${ genome_fasta.baseName }_and_ERCC92.fa"), emit: genome_fasta
        path("${ genome_gtf.baseName }_and_ERCC92.gtf"), emit: genome_gtf

    when:
        has_ercc

    script:
    """
    cat ${genome_fasta} ${ercc_fasta} > ${ genome_fasta.baseName }_and_ERCC92.fa
    cat ${genome_gtf} ${ercc_gtf} > ${ genome_gtf.baseName }_and_ERCC92.gtf
    """
}
