process DOWNLOAD_ERCC {
  storeDir "${reference_store_path}/ERCC_thermofisher"

  input:
    val(has_ercc)
    val(reference_store_path)

  output:
    path("ERCC92.fa") , emit: ercc_fasta
    path("ERCC92.gtf"), emit: ercc_gtf
    
  when:
    has_ercc

  script:
    """
    wget --no-check-certificate --quiet \
    -O ERCC92.zip \
    https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip  \
    && \
    unzip ERCC92.zip
    """
}