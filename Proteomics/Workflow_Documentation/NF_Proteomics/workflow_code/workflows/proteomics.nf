include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'


include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

ch_dp_tools_plugin = params.dp_tools_plugin ? 
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file("$projectDir/bin/dp_tools__NF_Proteomics"))

// Results go into output_dir/{OSD accession OR results_dir if provided}
// Default: run_dir/OSD-### or run_dir/results
// Create directory if it doesn't exist (checkIfExists: false)
output_dir = Channel.value(file(params.output_dir, type: 'dir', checkIfExists: false))



// Goal 1: get the output writing and the runsheet gen basic setup
workflow PROTEOMICS {
    take:
    main:
        // 
        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }
        
        // If accession is provided, parse accession input
        //   Use OSD accession to name the output directory unless params.results_dir is provided
        if ( params.accession ) {
            GET_ACCESSIONS( params.accession, params.api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            if (params.results_dir) {
                output_dir = output_dir.map { "$it/${params.results_dir}" }
            }
            else {
                output_dir = output_dir.combine(osd_accession).map { outdir, osd -> "$outdir/$osd" }
            }
        }
        else {
            if (params.results_dir) {
                output_dir = output_dir.map { "$it/${params.results_dir}" }
            }
            else {
                output_dir = output_dir.map { "$it/results" }
            }
        }
        output_dir = output_dir.first()

        // If no runsheet is provided, create one from the ISA archive or accession
        if ( params.runsheet == null ) {
            if ( params.isa_archive == null ) {
                FETCH_ISA( output_dir, osd_accession, glds_accession )
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, FETCH_ISA.out.isa_archive, ch_dp_tools_plugin )
                runsheet = ISA_TO_RUNSHEET.out.runsheet
            } else {
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, params.isa_archive, ch_dp_tools_plugin )
                runsheet = ISA_TO_RUNSHEET.out.runsheet
            }
        } else {
            runsheet = params.runsheet
        }

        // Convert the runsheet to a list of samples
        //samples = Channel.fromList(samplesheetToList(runsheet, "$projectDir/schema_runsheet.json"))

    emit:
        runsheet
}