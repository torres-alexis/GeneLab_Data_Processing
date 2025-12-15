// main.nf
nextflow.enable.dsl=2

// Command for: 'nextflow run main.nf --version'
if (params.version) {
    println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
    exit 0
}

include { PROTEOMICS } from './workflows/proteomics.nf'

// Main workflow
workflow {
    PROTEOMICS(
            )
}