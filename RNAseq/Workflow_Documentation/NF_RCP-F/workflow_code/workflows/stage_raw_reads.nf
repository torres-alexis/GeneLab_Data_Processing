nextflow.enable.dsl=2

workflow STAGE_RAW_READS {
    take:
        ch_samples
        stage_local
        truncate_to

    main:

        if ( stage_local && truncate_to ) {
        // Download truncated raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                 | branch {
                   paired: it.size() == 3
                   single: it.size() == 2
                 }
                 | set{ ch_raw_read_pointers }
            
            
            // TO DO: Move the two splitFastq truncation steps into processes that run in parallel. Low priority since this is just for debugging.

            // PAIRED END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                        // | view { it -> "TRUNCATED PAIRED READS ($truncate_to): $it[0]"}
                                        | set { ch_raw_reads }
            // SINGLE END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1] ] ]}
                                        // | view { it -> "TRUNCATED SINGLE READS ($truncate_to): $it[0]"}
                                        | mix( ch_raw_reads )
                                        | set { ch_raw_reads }

            // Moves the truncated files to expected raw read locations as per samplesheet
            ch_raw_reads | STAGE_READS
        } else if ( stage_local && !truncate_to ) {
        // download full raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                         | set { ch_raw_reads }

            // Download the raw reads and publish them to expected raw read locations as per samplesheet
            ch_raw_reads | STAGE_READS
        } else {
        // Don't download any raw reads
        }

    emit:
        raw_reads = stage_local ? STAGE_READS.out.raw_reads : null
}

process STAGE_READS {
    // Stages the raw reads into appropriate publish directory
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
