include { CLEAN_PATHS } from '../modules/clean_paths.nf'
include { PACKAGE_PROCESSING_INFO } from '../modules/package_processing_info.nf'
include { GENERATE_MD5SUMS } from '../modules/generate_md5sums.nf'
include { VALIDATE_PROCESSING } from '../modules/validate_processing.nf'

// Post-processing only: nextflow run main.nf --post_processing true ...
// Expected to only run after main workflow run.
// Expects processing_info/nextflow_processing_info_GLProteomics.txt, processing_info/nextflow_run_command_GLProteomics.txt, and processing_info/samples.txt in the processed output directory.
workflow POST_PROCESSING {
    main:
        processed_dir = "${params.output_dir}/${params.results_dir ?: (params.accession ?: 'results')}"
        ch_processed_directory = Channel.fromPath(processed_dir, type: 'dir', checkIfExists: true)
        ch_processing_info = Channel.fromPath("${processed_dir}/processing_info", type: 'dir', checkIfExists: true)

        if( params.clean_paths ) {
            CLEAN_PATHS( channel.value(processed_dir) )
            ch_processed_dir = CLEAN_PATHS.out.processed_dir
        } else {
            ch_processed_dir = channel.value(processed_dir)
        }

        ch_processing_info = ch_processing_info
            .combine( ch_processed_dir )
            .map { info, _dir -> info }
        ch_processed_directory = ch_processed_directory
            .combine( ch_processed_dir )
            .map { outdir, _dir -> outdir }

        PACKAGE_PROCESSING_INFO(ch_processing_info, processed_dir)
        GENERATE_MD5SUMS(ch_processed_directory, PACKAGE_PROCESSING_INFO.out.zip)
        VALIDATE_PROCESSING(
            ch_processed_directory,
            GENERATE_MD5SUMS.out.raw_md5sum,
            GENERATE_MD5SUMS.out.processed_md5sum
        )
}
