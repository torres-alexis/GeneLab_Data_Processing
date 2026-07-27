// main.nf

include { PROTEOMICS } from './workflows/proteomics.nf'
include { POST_PROCESSING as POST_PROCESSING_WORKFLOW } from './workflows/post_processing.nf'

workflow {
    if (params.version) {
        println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
        return
    }
    if (params.post_processing) {
        POST_PROCESSING_WORKFLOW()
    } else {
        PROTEOMICS()
    }
}
