include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { PARSE_RUNSHEET } from './parse_runsheet.nf'


// Set up channels
ch_dp_tools_plugin = params.dp_tools_plugin ? 
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file(params.mode == 'microbes' ? 
        "$projectDir/bin/dp_tools__NF_RCP_Bowtie2" : 
        "$projectDir/bin/dp_tools__NF_RCP"))
ch_glds = params.glds ? Channel.value(params.glds) : null
ch_osd = params.osd ? Channel.value(params.osd) : null
ch_runsheet = params.runsheet_path ? Channel.fromPath(params.runsheet_path) : null
ch_isa_archive = params.isa_archive_path ? Channel.fromPath(params.isa_archive_path) : null
ch_force_single_end = Channel.value(params.force_single_end)
ch_limit_samples_to = Channel.value(params.limit_samples_to)
workflow STAGE_ANALYSIS {
    // Takes in all required inputs to stage analysis
    // If runsheet path is provided, it will skip fetching ISA archive from OSDRq
    // Otherwise, it will fetch ISA archive from OSDR and convert to runsheet
    // It will then parse the runsheet to get sample metadata
    // It will then stage the raw reads - either full or truncated remote reads or local reads
    // It will then collect the sample IDs into a file
    // Emits:
    // - ch_all_raw_reads: all raw reads
    // - ch_samples_txt: New line separated list of sample IDs

    main:
        // If no runsheet path, fetch ISA archive from OSDR (if needed) and convert to runsheet
        if (ch_runsheet == null) {
            if (ch_isa_archive == null) {
                FETCH_ISA(ch_osd, ch_glds)
                ch_isa_archive = FETCH_ISA.out.isa_archive
            }
            ISA_TO_RUNSHEET( ch_glds, ch_isa_archive, ch_dp_tools_plugin )
            ch_runsheet = ISA_TO_RUNSHEET.out.runsheet
        }

        // Get sample metadata from runsheet
        PARSE_RUNSHEET( ch_runsheet, ch_force_single_end, ch_limit_samples_to )

        ch_samples = PARSE_RUNSHEET.out.samples

        // Stage full or truncated raw reads 
        STAGE_RAW_READS( ch_samples, params.stage_local, params.truncate_to )

        // Collect sample IDs into a file
        STAGE_RAW_READS.out.raw_reads | map{ it -> it[1] } | collect | set { ch_all_raw_reads }
        STAGE_RAW_READS.out.raw_reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_txt }

        // View the contents of ch_all_raw_reads
        // ch_all_raw_reads.view { it -> "All raw reads: $it" }

        // // View the contents of ch_samples_txt
        // ch_samples_txt.view { it -> "Samples txt file: $it" }

    emit:
        raw_reads = STAGE_RAW_READS.out.raw_reads
        ch_all_raw_reads // unused so far 
        samples_txt
}