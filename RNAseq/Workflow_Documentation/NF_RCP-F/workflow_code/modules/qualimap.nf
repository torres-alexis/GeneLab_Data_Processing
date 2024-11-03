process QUALIMAP_BAMQC {
    tag "Sample: ${ meta.id }"

    input:
    tuple val(meta), path(bam)
    path(feature_file) //gtf gff or bed
    val(strandedness)

    output:
    path("${meta.id}"), emit: results
    path  "versions.yml"       , emit: versions

    script:
    strandedness_opt_map = ["sense":"strand-specific-forward","antisense":"strand-specific-reverse","unstranded":"non-strand-specific"]
    def collect_pairs = meta.paired_end ? '--collect-overlap-pairs' : ''
    def memory = (task.memory.mega*0.8).intValue() + 'M'

    """
    unset DISPLAY
    mkdir -p tmp
    qualimap \\
        --java-mem-size=${memory} \\
        bamqc \\
        -bam ${bam} \\
        --feature-file ${feature_file} \\
        -p ${strandedness_opt_map.get(strandedness)} \\
        $collect_pairs \\
        -outdir ${meta.id} \\
        -nt ${task.cpus}

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    qualimap: \$(qualimap 2>&1 | sed 's/^.*QualiMap v.//; s/Built.*\$//')" >> versions.yml
    """
}