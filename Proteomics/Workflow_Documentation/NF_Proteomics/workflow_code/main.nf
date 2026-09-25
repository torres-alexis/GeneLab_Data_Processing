// main.nf

include { PROTEOMICS } from './workflows/proteomics.nf'
include { POST_PROCESSING as POST_PROCESSING_WORKFLOW } from './workflows/post_processing.nf'
include { validateParameters; paramsHelp } from 'plugin/nf-schema'

workflow {
    main:
        if (params.help) {
            println paramsHelp()
            exit 0
        }
        if (params.version) {
            println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
            ch_pub = Channel.empty()
        } else {
            validateParameters()
            if (params.post_processing) {
                POST_PROCESSING_WORKFLOW()
                ch_pub = POST_PROCESSING_WORKFLOW.out.published
            } else {
                PROTEOMICS()
                ch_pub = PROTEOMICS.out.published
            }
        }
    publish:
        published = ch_pub
}

output {
    published {
        path { dest, f -> f >> dest }
    }
}
