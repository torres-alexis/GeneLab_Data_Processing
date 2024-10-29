include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
// include { PREPARE_REFERENCES } from './prepare_references.nf'
include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

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
        reference_store_path
        derived_store_path
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
        
        has_ercc = ch_meta.map { it.has_ercc }

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
        gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.annotations_db_url

        // Use the parsed values if the reference_fasta and reference_gtf are not provided. Reference_source and reference_version are optional.
        if (reference_fasta == null && reference_gtf == null) {
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version

            DOWNLOAD_REFERENCES(reference_store_path, organism_sci, reference_fasta_url, reference_gtf_url, reference_version, reference_source)
            reference_fasta = DOWNLOAD_REFERENCES.out.reference_fasta_path // this is a path
            reference_gtf = DOWNLOAD_REFERENCES.out.reference_gtf_path // this is a path
        }
        else {
            reference_fasta.view {file -> "Using manually provided local reference genome fasta: ${file}"}
            reference_gtf.view {file -> "Using manually provided reference genome gtf: ${file}"}
        }
        
        // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
        if ( genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, reference_fasta, reference_gtf, reference_source, reference_version, genome_subsample )
            reference_fasta_pre_ercc = SUBSAMPLE_GENOME.out.subsampled_fasta
            reference_gtf_pre_ercc = SUBSAMPLE_GENOME.out.subsampled_gtf
        } else {
            reference_fasta_pre_ercc = reference_fasta
            reference_gtf_pre_ercc = reference_gtf
        }

        // Download and concatenate ERCC files if has_ercc is true
        if (has_ercc) {
            DOWNLOAD_ERCC(has_ercc, reference_store_path)
            ercc_fasta = DOWNLOAD_ERCC.out.ercc_fasta
            ercc_gtf = DOWNLOAD_ERCC.out.ercc_gtf
            CONCAT_ERCC(reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_pre_ercc, reference_gtf_pre_ercc, ercc_fasta, ercc_gtf, has_ercc)
            reference_fasta_post_ercc = CONCAT_ERCC.out.genome_fasta
            reference_gtf_post_ercc = CONCAT_ERCC.out.genome_gtf
        } else {
            ercc_fasta = null
            ercc_gtf = null
            reference_fasta_post_ercc = reference_fasta_pre_ercc
            reference_gtf_post_ercc = reference_gtf_pre_ercc
        }


        // Prepare reference files: reference genome with or without ERCC concatenation and annotation file
        // PREPARE_REFERENCES(reference_store_path, reference_fasta, reference_gtf, genome_subsample, has_ercc)
        // outputs:
        // - genome_fasta
        // - genome_gtf

    emit:
        reference_gtf_post_ercc
}
