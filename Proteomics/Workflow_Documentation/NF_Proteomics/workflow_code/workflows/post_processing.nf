include { PACKAGE_PROCESSING_INFO } from '../modules/package_processing_info.nf'
include { GENERATE_MD5SUMS } from '../modules/generate_md5sums.nf'
include { VALIDATE_PROCESSING } from '../modules/validate_processing.nf'
include { GENERATE_PROCESSED_PROTOCOL } from '../modules/generate_protocol.nf'

// Post-processing entry: nextflow run main.nf -entry POST_PROCESSING
// Expected to only run after main workflow run.
// Expects processing_scripts/nextflow_log_GLProteomics.txt, processing_scripts/nextflow_run_command_GLProteomics.txt, and processing_scripts/samples.txt in the processed output directory.
workflow POST_PROCESSING {
    main:
        processed_dir = "${params.output_dir}/${params.results_dir ?: (params.accession ?: 'results')}"
        ch_processed_directory = Channel.fromPath(processed_dir, type: 'dir', checkIfExists: true)
        ch_processing_info = Channel.fromPath("${processed_dir}/processing_scripts", type: 'dir', checkIfExists: true)
        ch_software_versions = Channel.fromPath("${processed_dir}/GeneLab/software_versions_*.md", checkIfExists: true)
        PACKAGE_PROCESSING_INFO(ch_processing_info, processed_dir)
        GENERATE_MD5SUMS(ch_processed_directory, PACKAGE_PROCESSING_INFO.out.zip)
        GENERATE_PROCESSED_PROTOCOL(ch_processed_directory, ch_software_versions)
        VALIDATE_PROCESSING(
            ch_processed_directory,
            GENERATE_MD5SUMS.out.raw_md5sum,
            GENERATE_MD5SUMS.out.processed_md5sum
        )
}
