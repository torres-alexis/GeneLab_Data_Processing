/*
* Workflow that accepts a GLDS accession and generates the following:
* 1. Download ISA.zip and generates RNASeq Samplesheet
* 2a. Downloads Raw Reads
* 2b. Downloads Truncated Raw Reads (Useful for testing with limited resources)
*/

// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

def mutate_to_single_end(it) {
  new_meta = it[0]
  new_meta["paired_end"] = false
  return [new_meta, it[1]]
}

// Import process from separate module file
include { GENERATE_METASHEET;
          STAGE_RAW_READS;
          get_runsheet_paths } from'./modules/genelab.nf'
include { RUNSHEET_FROM_GLDS } from './modules/RUNSHEET_FROM_GLDS.nf'
include { RUNSHEET_FROM_ISA } from './modules/RUNSHEET_FROM_ISA.nf'
/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow staging{
  take:
    ch_glds_accession
    stageLocal
  main:
    sample_limit = params.limitSamplesTo ? params.limitSamplesTo : -1 // -1 in take means no limit

    if (!params.runsheetPath && !params.isaArchivePath) {
      ch_osd_accession = Channel.from( params.osdAccession )
      RUNSHEET_FROM_GLDS(
        ch_osd_accession,
        ch_glds_accession,
        "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
      )
      RUNSHEET_FROM_GLDS.out.runsheet | set{ ch_runsheet }
      RUNSHEET_FROM_GLDS.out.isaArchive | set{ ch_isaArchive }
      GENERATE_METASHEET( RUNSHEET_FROM_GLDS.out.isaArchive, RUNSHEET_FROM_GLDS.out.runsheet )
    } else if (!params.runsheetPath && params.isaArchivePath) {
        ch_osd_accession = Channel.from( params.osdAccession )
        ch_isaArchive = Channel.from( params.isaArchivePath )
        RUNSHEET_FROM_ISA(
          ch_osd_accession,
          ch_glds_accession,
          ch_isaArchive,
          "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        )
        RUNSHEET_FROM_ISA.out.runsheet | set{ ch_runsheet }
        GENERATE_METASHEET( ch_isaArchive, RUNSHEET_FROM_ISA.out.runsheet )
    } else if ( params.runsheetPath && !params.isaArchivePath ) {
        ch_runsheet = channel.fromPath( params.runsheetPath )
    } else if ( params.runsheetPath && params.isaArchivePath ) {
        System.err.println("Error: User supplied both runsheetPath and isaArchivePath.  Only one or neither is allowed to be supplied!") // Print error message to System.err
        System.exit(1) // Exit with error code 1
    }

    ch_runsheet | splitCsv(header: true)
                | map{ row -> get_runsheet_paths(row) }
                | map{ it -> params.force_single_end ? mutate_to_single_end(it) : it }
                | take( sample_limit )
                | set{ ch_samples }

    if ( stageLocal && params.truncateTo ) {
      // download truncated raw reads
      // download full raw reads
      ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                 | branch {
                   paired: it.size() == 3
                   single: it.size() == 2
                 }
                 | set{ ch_raw_read_pointers}

       // PAIRED END
       // Only difference is the splitFastq arg 'pe'
      ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                  | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                  //| view { it -> "TRUNCATED PAIRED READS ($params.truncateTo): $it[0]"}
                                  | set { ch_raw_reads }
       // SINGLE END
       // Only difference is the splitFastq arg 'pe'
      ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: params.truncateTo, by: params.truncateTo, file: true)
                                  | map { it -> [ it[0], [ it[1] ] ]}
                                  //| view { it -> "TRUNCATED SINGLE READS ($params.truncateTo): $it[0]"}
                                  | mix( ch_raw_reads )
                                  | set { ch_raw_reads }

      // Moves the truncated files to expected raw read locations as per samplesheet
      ch_raw_reads | STAGE_RAW_READS

    } else if ( stageLocal && !params.truncateTo ) {
      // download full raw reads
      ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                 | set { ch_raw_reads }

      // Download the raw reads and publish them to expected raw read locations as per samplesheet
    ch_raw_reads | STAGE_RAW_READS

    } else {
      // Don't download any raw reads
    }


    emit:
      raw_reads = stageLocal ? STAGE_RAW_READS.out : null
      isa = params.isaArchivePath ? null : ch_isaArchive
      runsheet = ch_runsheet
      metasheet = params.runsheetPath ? null : GENERATE_METASHEET.out.metasheet
}