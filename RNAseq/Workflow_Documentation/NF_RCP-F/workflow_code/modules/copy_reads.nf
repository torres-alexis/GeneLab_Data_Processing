process COPY_READS {
    publishDir params.glds ? "${params.outdir}/GLDS-${params.glds}/00-RawData" : "${params.outdir}/00-RawData",
        mode: params.publish_dir_mode
    tag "${ meta.id }"

    input:
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}_R1_raw.fastq.gz
        cp -P 2.gz ${meta.id}_R2_raw.fastq.gz
        """
        } else {
        """
        cp -P 1.gz  ${meta.id}_raw.fastq.gz
        """
        }
}