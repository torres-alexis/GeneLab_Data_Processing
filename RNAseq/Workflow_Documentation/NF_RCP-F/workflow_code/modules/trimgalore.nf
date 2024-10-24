process TRIMGALORE {
    // Stages the raw reads into appropriate publish directory
    publishDir params.glds ? "${params.outdir}/GLDS-${params.glds}/01-TG_Preproc" : "${params.outdir}/01-TG_Preproc",
        mode: params.publish_dir_mode
  tag "${ meta.id }"
  label 'low_cpu_med_memory'

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
    path("versions.yml"), emit: versions

  script:

    """
    trim_galore --gzip \
    --cores $task.cpus \
    --phred33 \
    ${ meta.paired_end ? '--paired' : '' } \
    $reads \
    --output_dir .

    # rename with _trimmed suffix
    ${ meta.paired_end ? \
      "mv ${ meta.id }_R1_raw_val_1.fq.gz ${ meta.id }_R1_trimmed.fastq.gz; \
      mv ${ meta.id }_R2_raw_val_2.fq.gz ${ meta.id }_R2_trimmed.fastq.gz" : \
      "mv ${ meta.id }_raw_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz"}

    echo '"${task.process}":' > versions.yml
    echo "    trimgalore: \$(trim_galore -v | sed -n 's/.*version //p' | head -n 1)" >> versions.yml
    echo "    cutadapt: \$(cutadapt --version)" >> versions.yml
    """
}
