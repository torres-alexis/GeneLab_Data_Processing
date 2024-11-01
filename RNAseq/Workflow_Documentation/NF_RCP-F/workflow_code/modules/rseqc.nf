process ASSESS_STRANDEDNESS {

  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  stub:
    """
    assess_strandedness.py infer_out
    echo "unstranded:0.48595" > result.txt # override original results, this is because heavy truncation and genome subsampling can result in an ambigious strand assignment, which normally is an issue, but should be ignore for stubruns 
    """

  script:
    """
    assess_strandedness.py infer_out
    """
}

process GENEBODY_COVERAGE {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(genome_bed) // bam file sorted by coordinate

  output:
    path("${ meta.id }.geneBodyCoverage.txt"), emit: log_only
    path("${ meta.id }.geneBodyCoverage.*"), emit: all_output
    tuple val(meta), path("${ meta.id }.geneBodyCoverage.txt"), emit: log
    path("versions.txt"), emit: version

  script:
    def log_fname = "${ meta.id }.geneBodyCoverage.txt" 
    """    
    geneBody_coverage.py -r ${ genome_bed} -i ${ bam_file } -o ${ meta.id }

    # VERSIONS
    echo "RSeQC genebody_coverage version below:\n" > versions.txt 
    geneBody_coverage.py --version >> versions.txt
    """
}

process INFER_EXPERIMENT {
  tag "Sample: ${meta.id}"

  input:
    tuple val(meta), path(bam_file), path(bai_file) // bam file sorted by coordinate
    path(genome_bed)

  output:
    path("${meta.id}_infer_expt.out"), emit: log_only
    tuple val(meta), path("${meta.id}_infer_expt.out"), emit: log
    path("versions.yml"), emit: versions

  script:
    def log_fname = "${meta.id}_infer_expt.out"
    """    
    infer_experiment.py -r ${genome_bed} -i ${bam_file} -s ${params.rseqc_sample_count} > ${log_fname}

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    rseqc: \$(infer_experiment.py --version | sed -e 's/infer_experiment.py //g')" >> versions.yml
    """
}

process INNER_DISTANCE {
  tag "Sample: ${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(genome_bed) // bam file sorted by coordinate

  output:
    path("${ meta.id }.inner_distance_freq.txt"), emit: log_only
    path("${ meta.id }.inner_distance*"), emit: all_output
    tuple val(meta), path("${ meta.id }.inner_distance_freq.txt"), emit: log
    path("versions.txt"), emit: version

  when:
    meta.paired_end

  script:
    def log_fname = "${ meta.id }.inner_distance_freq.txt" 
    """    
    inner_distance.py -r ${ genome_bed } -i ${ bam_file } -k ${ params.quality.rseqc_sample_count } -l -150 -u 350 -o ${ meta.id } 

    # VERSIONS
    echo "RSeQC inner_distance version below:\n" > versions.txt 
    inner_distance.py --version >> versions.txt
    """
}

process READ_DISTRIBUTION {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(genome_bed) // bam file sorted by coordinate

  output:
    path("${ meta.id }_read_dist.out"), emit: log_only
    tuple val(meta), path("${ meta.id }_read_dist.out"), emit: log
    path("versions.txt"), emit: version

  script:
    def log_fname = "${ meta.id }_read_dist.out"
    """    
    read_distribution.py -r ${ genome_bed } -i ${ bam_file } > ${ meta.id }_read_dist.out

    # VERSIONS
    echo "RSeQC read_distribution version below:\n" > versions.txt 
    read_distribution.py --version >> versions.txt
    """
}