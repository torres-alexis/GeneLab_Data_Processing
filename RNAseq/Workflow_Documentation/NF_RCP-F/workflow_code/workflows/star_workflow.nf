include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from './parse_annotations_table.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
workflow STAR_WORKFLOW {
    take:
        dp_tools_plugin
        annotations_csv_url_string
        force_single_end
        limit_samples_to
        accession
        isa_archive_path
        runsheet_path
        api_url
    main:
        // Set up runsheet
        if (runsheet_path == null) {
            //Get both OSD and GLDS accessions based on the input accession
            GET_ACCESSIONS( accession, api_url )
            accessions_txt = GET_ACCESSIONS.out.accessions_txt // returns accessions.txt with line1 = osd_accession, line2 = glds_accession. 
            osd_accession = accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = accessions_txt.map { it.readLines()[1].trim() }

            //Fetch ISA archive if not provided
            if (isa_archive_path == null) {
                FETCH_ISA(osd_accession, glds_accession)
                isa_archive = FETCH_ISA.out.isa_archive
            }
            //Convert ISA archive to runsheet
            ISA_TO_RUNSHEET( osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
        }

        //Parse runsheet
        PARSE_RUNSHEET(runsheet_path)

    emit:
        accessions_txt
}
