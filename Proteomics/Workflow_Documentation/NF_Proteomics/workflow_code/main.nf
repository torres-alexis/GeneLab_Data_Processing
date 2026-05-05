// main.nf

// Command for: 'nextflow run main.nf --version'
if (params.version) {
    println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
    exit 0
}

include { PROTEOMICS } from './workflows/proteomics.nf'
include { POST_PROCESSING as POST_PROCESSING_WORKFLOW } from './workflows/post_processing.nf'

// Main workflow (default entry)
workflow {
    PROTEOMICS()
}

// Post-processing entry: nextflow run main.nf -entry POST_PROCESSING
workflow POST_PROCESSING {
    POST_PROCESSING_WORKFLOW()
}
