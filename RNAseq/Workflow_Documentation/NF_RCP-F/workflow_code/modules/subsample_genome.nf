process SUBSAMPLE_GENOME {
  tag "Subsample: Region ${genome_subsample}"
  storeDir "${derived_store_path}/subsampled_files/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"
  input:
    val(derived_store_path)
    val(organism_sci)
    path(reference_fasta)
    path(reference_gtf)
    val(reference_source)
    val(reference_version)
    val(genome_subsample)

  output:
    path("${ reference_fasta.baseName }_sub_${ genome_subsample  }.fa"), emit: subsampled_fasta
    path("${ reference_gtf.baseName }_sub_${ genome_subsample }.gtf"), emit: subsampled_gtf

  script:
    """
    # Extract sequence using output flag instead of redirection
    samtools faidx ${reference_fasta} ${genome_subsample} -o ${ reference_fasta.baseName }_sub_${ genome_subsample }.fa

    # subsample gtf file
    grep '^#!' ${reference_gtf } > ${ reference_gtf.baseName }_sub_${ genome_subsample  }.gtf
    grep '^${genome_subsample}\t' ${reference_gtf} >> ${ reference_gtf.baseName }_sub_${ genome_subsample  }.gtf
    """
}