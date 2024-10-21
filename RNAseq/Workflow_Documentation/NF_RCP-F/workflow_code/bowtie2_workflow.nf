// color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { FASTQC as RAW_FASTQC } from './modules/quality.nf'
include { FASTQC as TRIMMED_FASTQC } from './modules/quality.nf'
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"raw")
include { MULTIQC as TRIMMED_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"trimmed")
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"trimming")
include { MULTIQC as ALIGN_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"align")
def countsLabel = params.microbes ? "FeatureCounts" : "RSEM_count"
include { MULTIQC as COUNT_MULTIQC } from './modules/quality.nf' addParams(MQCLabel: countsLabel)
include { MULTIQC as ALL_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"all")
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_BOWTIE2;
          ALIGN_BOWTIE2;
          SUBSAMPLE_GENOME;
          CONCAT_ERCC;
          QUANTIFY_BOWTIE2_GENES;
          QUANTIFY_NONZERO_GENES } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/DGE_BY_DESEQ2'
include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_BOWTIE2_ALIGNMENTS;
          VV_RSEQC;
          VV_FEATURECOUNTS;
          VV_DESEQ2_ANALYSIS;
          VV_CONCAT_FILTER } from './modules/vv.nf'
include { GET_MAX_READ_LENGTH } from './modules/fastqc.nf'
include { UPDATE_ISA_TABLES;
          GENERATE_MD5SUMS;
          SOFTWARE_VERSIONS } from './modules/genelab.nf'

include { staging as STAGING } from './stage_analysis.nf'
include { references as REFERENCES } from './references.nf'
include { strandedness as STRANDEDNESS } from './strandedness.nf'

// Set up channel containing glds accession number
if ( params.gldsAccession ) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}
// Set up channel containing osd accession number
if ( params.osdAccession ) {ch_osd_accession = Channel.from( params.osdAccession )}
// Set the MultiQC config channel
ch_multiqc_config = params.multiqcConfig ? Channel.fromPath(params.multiqcConfig) : Channel.fromPath("NO_FILE")

workflow BOWTIE2_WORKFLOW  {
	main:
    STAGING( ch_glds_accession, params.stageLocal )
    // This process can use a single meta and a collection of read paths
    STAGING.out.raw_reads | first 
                          | map{it -> it[0]} 
                          | view { meta -> "${c_bright_green}Autodetected Processing Metadata:\n\t hasERCC: ${meta.has_ercc}\n\t pairedEND: ${meta.paired_end}\n\t organism: ${meta.organism_sci}${c_reset}"  }
                          | set { ch_meta }

    STAGING.out.raw_reads | map{ it -> it[1] } | collect | set { ch_all_raw_reads }
    STAGING.out.raw_reads | map { it[0].id }
                          | collectFile(name: "samples.txt", sort: true, newLine: true)
                          | set { ch_samples_txt }

    STAGING.out.raw_reads | RAW_FASTQC

    RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
                          | flatten
                          | unique
                          | collect
                          | set { raw_mqc_ch }
    RAW_FASTQC.out.fastqc | map { it -> [ it[2] ] }
                          | flatten
                          | GET_MAX_READ_LENGTH

    GET_MAX_READ_LENGTH.out.length  | max { it.toInteger() }
                                    | set { max_read_length_ch }

    STAGING.out.raw_reads |  TRIMGALORE

    TRIMGALORE.out.reads | TRIMMED_FASTQC

    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } \
                              | flatten \
                              | unique \
                              | collect \
                              | set { trim_mqc_ch }

    REFERENCES( ch_meta | map { it.organism_sci }, ch_meta | map { it.has_ercc } )
    REFERENCES.out.genome_annotations | set { genome_annotations }      

    BUILD_BOWTIE2(
      genome_annotations,
      ch_meta,
      REFERENCES.out.reference_version_and_source
    )
    TRIMGALORE.out.reads | combine( BUILD_BOWTIE2.out.build ) | ALIGN_BOWTIE2

    ALIGN_BOWTIE2.out.alignment_logs | collect | set { align_mqc_ch }

    STRANDEDNESS ( ALIGN_BOWTIE2.out.bam, REFERENCES.out.genome_bed, ch_samples_txt )
    STRANDEDNESS.out.strandedness | map { it.text.split(":")[0] } | set { strandedness_ch }

    bam_paths = ALIGN_BOWTIE2.out.bam.map { it[1] } | toSortedList()

    QUANTIFY_BOWTIE2_GENES(
      ch_meta,
      genome_annotations,
      ch_samples_txt,
      strandedness_ch,
      bam_paths
    )
      
    QUANTIFY_BOWTIE2_GENES.out.publishables
      .map { it[0] }
      .set { bowtie2_gene_counts }

    QUANTIFY_BOWTIE2_GENES.out.publishables
      .map { it[1] }
      .set { fcsummary_ch }

    QUANTIFY_NONZERO_GENES(bowtie2_gene_counts) | set {numNonzeroGenes}


    // Note: This is loaded as a zip file to ensure correct caching (directories don't seem to yield identical hases)
    ch_r_scripts = channel.fromPath( "${ projectDir }/bin/dge_annotation_R_scripts.zip" ) 
  
    DGE_BY_DESEQ2(STAGING.out.runsheet, 
                  bowtie2_gene_counts, 
                  ch_meta, 
                  REFERENCES.out.gene_annotations, 
                  ch_r_scripts
                  )

    // ALL MULTIQC
    RAW_MULTIQC( ch_samples_txt, raw_mqc_ch, ch_multiqc_config  )
    TRIMMED_MULTIQC( ch_samples_txt, trim_mqc_ch, ch_multiqc_config ) // refering to the trimmed reads
    TRIM_MULTIQC( ch_samples_txt, TRIMGALORE.out.reports | collect, ch_multiqc_config ) // refering to the trimming process
    ALIGN_MULTIQC( ch_samples_txt, align_mqc_ch, ch_multiqc_config )
    COUNT_MULTIQC( ch_samples_txt, fcsummary_ch, ch_multiqc_config )
      raw_mqc_ch | concat( trim_mqc_ch ) 
        | concat( ALIGN_BOWTIE2.out.alignment_logs ) 
        | concat( STRANDEDNESS.out.rseqc_logs )
        | concat( fcsummary_ch )
        | concat( TRIMGALORE.out.reports )
        | collect | set { all_mqc_ch }
    ALL_MULTIQC( ch_samples_txt, all_mqc_ch, ch_multiqc_config )

    VV_RAW_READS( ch_meta,
                  STAGING.out.runsheet,
                  ch_all_raw_reads,
                  RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } | flatten | collect,
                  RAW_MULTIQC.out.zipped_report,
                  RAW_MULTIQC.out.unzipped_report,
                  "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
                )
    VV_TRIMMED_READS( ch_meta,
                      STAGING.out.runsheet,
                      TRIMGALORE.out.reads | map { it -> it[1] } | flatten | collect,
                      TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } | flatten | collect,
                      TRIMMED_MULTIQC.out.zipped_report,
                      TRIMMED_MULTIQC.out.unzipped_report,
                      TRIMGALORE.out.reports | collect,
                      TRIM_MULTIQC.out.zipped_report,
                      TRIM_MULTIQC.out.unzipped_report,
                      "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
                    )
    VV_BOWTIE2_ALIGNMENTS( STAGING.out.runsheet,
                      ALIGN_BOWTIE2.out.publishables | collect,
                      ALIGN_MULTIQC.out.zipped_report,
                      ALIGN_MULTIQC.out.unzipped_report,
                      STRANDEDNESS.out.bam_bed | collect,
                      "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
    )
    VV_FEATURECOUNTS( STAGING.out.runsheet,
                        QUANTIFY_BOWTIE2_GENES.out.publishables | collect,
                        numNonzeroGenes,
                        COUNT_MULTIQC.out.zipped_report,
                        COUNT_MULTIQC.out.unzipped_report,
                        "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
    )
    VV_RSEQC( ch_meta,
                  STAGING.out.runsheet,
                  STRANDEDNESS.out.rseqc_logs,
                  STRANDEDNESS.out.genebody_coverage_multiqc,
                  STRANDEDNESS.out.infer_experiment_multiqc,
                  STRANDEDNESS.out.inner_distance_multiqc,
                  STRANDEDNESS.out.read_distribution_multiqc,
                  "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
    )

    VV_DESEQ2_ANALYSIS( ch_meta,
                        STAGING.out.runsheet,
                        QUANTIFY_BOWTIE2_GENES.out.publishables,
                        COUNT_MULTIQC.out.zipped_report,
                        COUNT_MULTIQC.out.unzipped_report,
                        DGE_BY_DESEQ2.out.norm_counts,
                        DGE_BY_DESEQ2.out.dge,
                        DGE_BY_DESEQ2.out.norm_counts_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
                        DGE_BY_DESEQ2.out.dge_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
                        "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" // dp_tools plugin
    )

    // Software Version Capturing
    nf_version = "Nextflow Version:".concat("${nextflow.version}\n<><><>\n")
    ch_nextflow_version = Channel.value(nf_version)
    ch_software_versions = Channel.empty()
    RAW_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    RAW_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMGALORE.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMMED_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMMED_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    ALIGN_BOWTIE2.out.version | mix(ch_software_versions) | set{ch_software_versions}
    QUANTIFY_BOWTIE2_GENES.out.version | mix(ch_software_versions) | set{ch_software_versions}
    DGE_BY_DESEQ2.out.version | mix(ch_software_versions) | set{ch_software_versions}
    STRANDEDNESS.out.versions | mix(ch_software_versions) | set{ch_software_versions}
    ch_software_versions | map { it.text + "\n<><><>\n"}
                          | unique
                          | mix(ch_nextflow_version)
                          | collectFile(name: "software_versions.txt", newLine: true, cache: false)
        | set{ch_final_software_versions}

    // VV processes
      
    VV_CONCAT_FILTER( VV_RAW_READS.out.log | mix( VV_TRIMMED_READS.out.log,
                                                  VV_BOWTIE2_ALIGNMENTS.out.log,
                                                  VV_RSEQC.out.log,
                                                  VV_FEATURECOUNTS.out.log,
                                                  VV_DESEQ2_ANALYSIS.out.log,
                                                  ) | collect )

    // Generate final versions output
    SOFTWARE_VERSIONS(ch_final_software_versions)

    emit:
      meta = ch_meta 
      runsheet = STAGING.out.runsheet 
      done = SOFTWARE_VERSIONS.out 
}