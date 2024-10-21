nextflow.enable.dsl=2

// color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { staging as STAGING } from './stage_analysis.nf'
// Include post-processing modules
include { UPDATE_ISA_TABLES;
          GENERATE_MD5SUMS;
          SOFTWARE_VERSIONS } from './modules/genelab.nf'

// Conditionally include star or bowtie workflow
if( params.microbes ) {
    include { BOWTIE2_WORKFLOW } from './bowtie2_workflow.nf'
} else {
    include { STAR_WORKFLOW } from './star_workflow.nf'
}

// Help menu
if (params.help) {
    println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
    println("┇ RNASeq Consensus Pipeline: $workflow.manifest.version  ┇")
    println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
    println("Usage example 1: Processing GLDS datasets using genome fasta and gtf from Ensembl")
    println("   > nextflow run ./main.nf --gldsAccession GLDS-194")
    println()
    println("Usage example 2: Processing GLDS datasets using local genome fasta and gtf")
    println("   Note: ensemblVersion and ref_source are used here to label subdirectories for derived reference files.")
    println("   > nextflow run ./main.nf --osdAccession OSD-194 --gldsAccession GLDS-194 --ensemblVersion 96 --ref_source <reference_label>  --ref_fasta </path/to/fasta> --ref_gtf </path/to/gtf>")
    println()
    println("Usage example 3: Processing Other datasets")
    println("   Note: This requires a user-created runsheet.")
    println("   > nextflow run ./main.nf --runsheetPath </path/to/runsheet>")
    println()
    println("arguments:")
    println("  --help                show this help message and exit")
    println(" --microbes             Use this argument to use Bowtie2 and FeatureCounts to align and quantify genes instead of using STAR and RSEM.")
    println("  --osdAccession OSD-000")
    println("                        the OSD accession id to process through the RNASeq Consensus Pipeline.")
    println("  --gldsAccession GLDS-000")
    println("                        the GLDS accession id to process through the RNASeq Consensus Pipeline.")
    println("  --runsheetPath        Use a local runsheet instead one automatically generated from a GLDS ISA archive.")
    println("  --ensemblVersion n    Specifies the ensembl Version to use for the reference genome. The default version is ")
    println("  --skipVV              Skip automated V&V. Default: false")
    println("  --outputDir           Directory to save staged raw files and processed files. Default: <launch directory>")
    println("  --limitSamplesTo n    limit the number of samples staged to a number.")
    println("  --genomeSubsample n   subsamples genome fasta and gtf files to the supplied chromosome.")
    println("  --truncateTo n        limit number of reads downloaded and processed to *n* reads , for paired end limits number of reverse and forward read files to *n* reads each.")
    println("  --force_single_end    forces analysis to use single end processing.  For paired end datasets, this means only R1 is used.  For single end studies, this should have no effect.")
    println("  --stageLocal          download the raw reads files for the supplied GLDS accession id.  Set to false to retrieve metadata and generate a runsheet for GLDS datasets to disable raw read download and processing.  Default: true")
    println("  --ref_fasta           specifies a reference fasta from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl")  
    println("  --ref_gtf             specifies a reference gtf from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl")  
    println("  --referenceStorePath  specifies the directory where fetched reference files are downloaded to")  
    println("  --derivedStorePath    specifies the directory where derivative reference files are saved. Examples of such files in this pipeline included BED and PRED files generated from the reference gtf")  
    println("  --ref_source          a string to label subdirectories in 'StorePath' paths. Examples include 'ensembl' or 'ensembl_plants'.")  
    println("  -stub-run             runs the workflow forcing 'unstranded' RSEM settings and using dummy gene counts in the differential gene expression (DGE) analysis. Useful when combined with the --truncateTo parameter this often leads to low gene counts and errors in the DGE analysis")  
    exit 0
}

println "PARAMS: $params"
println "\n"
println "Storing any newly fetched primary references files here: ${params.referenceStorePath}"
println "Storing any newly generated derived reference files here: ${params.derivedStorePath}"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( params.gldsAccession ) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}

// Check conditionally required parameter (if using direct fasta, an ensemblVersion must also be supplied)
if ( params.ref_fasta ) {
  if ( !params.ensemblVersion ) { exit 1, "Missing Required Parameter: ensemblVersion. Example for setting on CLI: --ensemblVersion 96" }
}

if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")

/**************************************************
* DEBUG WARNING  **********************************
**************************************************/
if ( params.limitSamplesTo || params.truncateTo || params.force_single_end || params.genomeSubsample) {
  println("${c_back_bright_red}WARNING WARNING: DEBUG OPTIONS ENABLED!")
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.genomeSubsample ? println("Subsampling reference genome to chromosome '${params.genomeSubsample}'") : println("No subsampling of reference genome")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
  println("WARNING WARNING: DEBUG OPTIONS ENABLED!${c_reset}")
} else {
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.genomeSubsample ? println("Subsampling reference genome to chromosome '${params.genomeSubsample}'") : println("No subsampling of reference genome")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
}

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/
if ( params.stageLocal && params.truncateTo ) {
  // download truncated raw reads
  println("${c_bright_green}Staging truncated raw reads for ${params.gldsAccession}${c_reset}")
} else if ( params.stageLocal && !params.truncateTo ) {
  // download full raw reads
  println("${c_bright_green}Staging raw reads for ${params.gldsAccession}${c_reset}")
} else {
  // maybe print some nice data from the samplesheet
  println("${c_bright_green}No Staging of raw reads.  Only getting Metadata for ${params.gldsAccession}${c_reset}")
}

// Workflow selection based on params
workflow {
    if (params.microbes) {
        BOWTIE2_WORKFLOW()
    } else {
        STAR_WORKFLOW()
    }
}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.outputDir }/${ params.gldsAccession }"
      println "V&V logs location: ${ params.outputDir }/${ params.gldsAccession }/VV_Logs"
      println "Pipeline tracing/visualization files location:  ${ params.outputDir }/${ params.gldsAccession }/Resource_Usage${c_reset}"
    }
}

workflow STAGING_ONLY {
  main:
    STAGING( ch_glds_accession, false )
}

workflow POST_PROCESSING {
  main:
    ch_processed_directory = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }", checkIfExists: true)
    ch_runsheet = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }/Metadata/*_runsheet.csv", checkIfExists: true)
    // Conditionally set the dp_tools directory based on the aligner method
    ch_dptools_dir = Channel.fromPath(params.microbes ? "${ projectDir }/bin/dp_tools__NF_RCP_Bowtie2" : "${ projectDir }/bin/dp_tools__NF_RCP")
    GENERATE_MD5SUMS(ch_processed_directory, ch_runsheet, ch_dptools_dir)
    UPDATE_ISA_TABLES(ch_processed_directory, ch_runsheet, ch_dptools_dir)
}