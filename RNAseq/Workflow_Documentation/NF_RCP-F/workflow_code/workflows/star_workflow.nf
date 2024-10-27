include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
workflow STAR_WORKFLOW {
    take:
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        force_single_end
        limit_samples_to
        truncate_to
        genome_subsample
        reference_source
        reference_version
        reference_fasta
        reference_gtf
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

        //Parse runsheet for processing metadata and check for unique read file paths
        PARSE_RUNSHEET(runsheet_path)
        // samples is a channel of tuples, where each tuple contains:
        // 1. A map of metadata about the sample with the following keys:
        //    - id: String (Sample Name)
        //    - organism_sci: String (Organism name, lowercase with underscores, e.g. 'mus_musculus')
        //    - paired_end: Boolean
        //    - has_ercc: Boolean
        //    - factors: Map of factor values (keys are factor names, values are factor levels)
        // 2. A list of read file paths associated with that sample
        // Example usage: samples.take(1) | view { meta, reads -> ... }
        samples = PARSE_RUNSHEET.out.samples

        // Extract metadata from the first sample and set it as a channel
        samples | first 
                | map { meta, reads -> meta }
                | set { ch_meta }
        
        // // Add a view operation to display the metadata
        // ch_meta | view { meta -> 
        //     """
        //     Metadata for first sample:
        //     Sample ID: ${meta.id}
        //     Organism: ${meta.organism_sci}
        //     Paired End: ${meta.paired_end}
        //     Has ERCC: ${meta.has_ercc}
        //     Factors: ${meta.factors}
        //     """
        // }



        // Extract organism_sci from the first sample in order to read the correct row from the annotations table
        ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }

        // Add a view operation to display the organism_sci
        // organism_sci | view { it -> 
        //     "Organism: ${it}"
        // }


        // Parse the annotations table
        // Outputs:
        // - reference_fasta_url: val(reference_fasta_url)
        // - reference_gtf_url: val(reference_gtf_url)
        // - annotations_db_url: val(annotations_db_url)
        // - reference_source: val(reference_source)
        // - reference_version: val(reference_version)
        PARSE_ANNOTATIONS_TABLE(annotations_csv_url_string, organism_sci)
        // Use the parsed values if the reference_fasta and reference_gtf are not provided. Reference_source and reference_version are optional.
        if (reference_fasta == null || reference_gtf == null) {
            reference_fasta = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
        }
        
        // Add a view operation to display the reference_fasta
        // reference_fasta.view { it -> println "Reference FASTA URL: ${it}" }

    emit:
        accessions_txt
}
