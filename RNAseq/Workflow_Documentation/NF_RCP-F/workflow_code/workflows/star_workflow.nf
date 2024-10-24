include { STAGE_ANALYSIS } from './stage_analysis.nf'
include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { MULTIQC as RAW_MULTIQC } from '../modules/multiqc.nf'  addParams(MQCLabel:"raw")
include { MULTIQC as TRIMMED_MULTIQC } from '../modules/multiqc.nf'  addParams(MQCLabel:"trimmed")
workflow STAR_WORKFLOW {
    // Staging the analysis and set up input channels
    STAGE_ANALYSIS()
    ch_raw_reads = STAGE_ANALYSIS.out.raw_reads
    ch_samples_txt = STAGE_ANALYSIS.out.samples_txt

    RAW_FASTQC(ch_raw_reads)

    RAW_FASTQC.out.zip
        .map { meta, zip -> zip }  // Keep only the fastqc output zip files 
        .collect()                 // Collect files into a list
        .set { ch_raw_fastqc_zip } 

    ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")
    RAW_MULTIQC( ch_samples_txt, ch_raw_fastqc_zip, ch_multiqc_config  )

    // Get the maximum read length
    GET_MAX_READ_LENGTH(ch_raw_fastqc_zip)
    max_read_length = GET_MAX_READ_LENGTH.out.length

    // View to see the maximum read length
    // max_read_length.view { "Max read length: $it" }

    STAGE_ANALYSIS.out.raw_reads |  TRIMGALORE

    TRIMGALORE.out.reads | TRIMMED_FASTQC
    TRIMMED_FASTQC.out.zip
        .map { meta, zip -> zip }          // Keep only the fastqc output zip files 
        .collect()                         // Collect files into a list
        .set { ch_trimmed_fastqc_zip }    

    
}