process DOWNLOAD_REFERENCES {
  // Download and decompress genome and gtf files
  tag "Organism: ${organism_sci} | Reference Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ' | Reference Version: ' + reference_version : ''}"
  storeDir "${reference_store_path}/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"

  input:
    val(reference_store_path)
    val(organism_sci)
    val(fasta_url)
    val(gtf_url)
    val(reference_version)
    val(reference_source)
  
  output:
    path("{*.fa,*.fna}"), emit: reference_fasta_path //
    path("*.gtf")       , emit: reference_gtf_path

  script:
  """
  wget ${fasta_url}
  gunzip *.gz

  wget ${gtf_url}
  gunzip *.gz
  """
}