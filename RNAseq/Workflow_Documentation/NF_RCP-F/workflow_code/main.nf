// main.nf

nextflow.enable.dsl=2

def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

// Include the showHelp function from help.nf
include { showHelp } from './modules/help.nf'
// Only display the help message if --help parameter is specified
if (params.help) {
    showHelp(workflow)
}

// Print the pipeline version
println """
${colorCodes.c_bright_green}${colorCodes.c_line}
GeneLab RNA-Seq Consensus Pipeline NF-RCP-${workflow.manifest.version}
${colorCodes.c_line}${colorCodes.c_reset}
""".stripIndent()

// Debug warning
println("${colorCodes.c_yellow}")
if (params.limit_samples_to || params.truncate_to || params.force_single_end || params.genome_subsample) {
    println("WARNING: Debugging options enabled!")
    println("Sample limit: ${params.limit_samples_to ?: 'Not set'}")
    println("Read truncation: ${params.truncate_to ? "First ${params.truncate_to} records" : 'Not set'}")
    println("Reference genome subsampling: ${params.genome_subsample ? "Chromosome '${params.genome_subsample}'" : 'Not set'}")
    println("Force single-end analysis: ${params.force_single_end ? 'Yes' : 'No'}")
} else {
    println("No debugging options enabled")
}
println("${colorCodes.c_reset}")

// Check required parameters
if ((params.glds && params.osd) || params.runsheet_path || params.isa_archive_path) {
    // Proceed
} else {
    log.error """
        Missing Required Parameters: You must provide either both --osd and --glds, or --runsheet_path, or --glds and --isa_archive_path.
        Examples:
          --osd 194 --glds 194
          --runsheet_path /path/to/runsheet.csv
          --isa_archive_path /path/to/isa_archive.zip
    """
    exit 0
}

include { STAR_WORKFLOW } from './workflows/star_workflow.nf'

// Main workflow
workflow {
    if (params.mode == 'microbes') {
        //BOWTIE2_WORKFLOW()
    } else {
        STAR_WORKFLOW()
    }
}
