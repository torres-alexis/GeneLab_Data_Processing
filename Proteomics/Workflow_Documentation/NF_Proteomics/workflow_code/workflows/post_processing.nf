include { CLEAN_PATHS } from '../modules/clean_paths.nf'
include { PACKAGE_PROCESSING_INFO } from '../modules/package_processing_info.nf'
include { GENERATE_MD5SUMS } from '../modules/generate_md5sums.nf'
include { VALIDATE_PROCESSING } from '../modules/validate_processing.nf'

def as_list(x) {
    if (x == null) return []
    if (x instanceof Collection && !(x instanceof CharSequence)) return x as List
    return [x]
}

def dest_file(root, dir, f) {
    return ["${root}/${dir}/${f.name}".toString(), f]
}

def pub(ch, root, dir) {
    return ch.combine(root).combine(channel.value(dir)).flatMap { row ->
        def items = as_list(row)
        def d = items[-1]
        def r = items[-2]
        as_list(items[0]).collect { f -> dest_file(r, d, f) }
    }
}

// nextflow run main.nf --post_processing true ...
// Needs processing_info/{nextflow_processing_info,nextflow_run_command,samples}* after the main run.
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

        ch_root = Channel.value((params.results_dir ?: (params.accession ?: 'results')).toString())
        ch_published = pub(PACKAGE_PROCESSING_INFO.out.zip, ch_root, 'GeneLab')
            .mix(pub(GENERATE_MD5SUMS.out.raw_md5sum, ch_root, 'GeneLab'))
            .mix(pub(GENERATE_MD5SUMS.out.processed_md5sum, ch_root, 'GeneLab'))
            .mix(pub(VALIDATE_PROCESSING.out.validation_log, ch_root, 'GeneLab'))

    emit:
        published = ch_published
}
