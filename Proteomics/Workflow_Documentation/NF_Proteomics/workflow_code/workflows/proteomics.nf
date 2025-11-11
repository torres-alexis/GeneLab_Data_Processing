include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { STAGE_INPUT } from '../modules/stage_input.nf'
include { RAWBEANS_QC } from '../modules/rawbeans_qc.nf'
include { RAWBEANS_QC_ALL } from '../modules/rawbeans_qc.nf'
include { OPENMS_FILEINFO } from '../modules/openms_fileinfo.nf'
include { COMET_SEARCH } from '../modules/comet_search.nf'
include { MSFRAGGER } from '../modules/msfragger.nf'
include { MSFRAGGER_CONFIG } from '../modules/msfragger.nf'
include { GET_PROTEOME } from '../modules/get_proteome.nf'
include { CRYSTAL_C } from '../modules/crystal_c.nf'
include { CRYSTAL_C_CONFIG } from '../modules/crystal_c.nf'
include { PHILOSOPHER } from '../modules/philosopher.nf'
include { IONQUANT } from '../modules/ionquant.nf'
include { MSSTATS } from '../modules/msstats.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

ch_dp_tools_plugin = params.dp_tools_plugin ? 
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file("$projectDir/bin/dp_tools__NF_Proteomics"))

// Results go into output_dir/{OSD accession OR results_dir if provided}
// Default: run_dir/OSD-### or run_dir/results
// Create directory if it doesn't exist (checkIfExists: false)
output_dir = Channel.value(file(params.output_dir, type: 'dir', checkIfExists: false))



// Goal 1: get the output writing and the runsheet gen basic setup
workflow PROTEOMICS {
    take:
    main:
        // 
        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }
        ch_reference_store_path = Channel.value(params.reference_store_path)
        
        // Handle accession and output directory naming
        if ( params.accession ) {
            GET_ACCESSIONS( params.accession, params.api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            
            // Use OSD accession for output directory unless params.results_dir is provided
            if (params.results_dir) {
                output_dir = output_dir.map { "$it/${params.results_dir}" }
            } else {
                output_dir = output_dir.combine(osd_accession).map { outdir, osd -> "$outdir/$osd" }
            }
        } else {
            // No accession provided, use results_dir or default
            if (params.results_dir) {
                output_dir = output_dir.map { "$it/${params.results_dir}" }
            } else {
                output_dir = output_dir.map { "$it/results" }
            }
        }

        // Handle runsheet: generate from ISA or use provided
        if ( params.runsheet_path == null ) {
            // Generate runsheet from ISA archive
            if ( params.isa_archive == null ) {
                FETCH_ISA( output_dir, osd_accession, glds_accession )
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, FETCH_ISA.out.isa_archive, ch_dp_tools_plugin )
                runsheet = ISA_TO_RUNSHEET.out.runsheet
            } else {
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, params.isa_archive, ch_dp_tools_plugin )
                runsheet = ISA_TO_RUNSHEET.out.runsheet
            }
        } else {
            // Use provided runsheet
            runsheet = params.runsheet_path
        }

        // Convert the runsheet to a list of samples and run processing
        samples = Channel.fromList(samplesheetToList(runsheet, "$projectDir/schema_runsheet.json"))
            .map { row -> 
                def meta = [
                    id: row[0],           // Sample Name
                ]
                return tuple(meta, row[1])  // meta, input_file
            }
        samples.view { meta, input_file ->
            """Sample: ${meta.id}, Input File: ${input_file}"""
        }

        // Stage input files for each sample
        STAGE_INPUT(output_dir, samples)
        // Run QC on each sample's raw data
        //RAWBEANS_QC(output_dir, STAGE_INPUT.out.mzml_files)
        // Run QC on all samples' raw data
        RAWBEANS_QC_ALL(output_dir, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())

        // Validate database input: cannot specify both uniprot_id and reference_proteome
        if (params.uniprot_id && params.reference_proteome) {
            error "ERROR: Cannot specify both uniprot_id and reference_proteome. Use EITHER uniprot_id to download from UniProt OR reference_proteome to use a custom FASTA file."
        }
        
        if (!params.uniprot_id && !params.reference_proteome) {
            error "ERROR: Must specify either uniprot_id or reference_proteome."
        }
        
        // Get proteome database: download from UniProt OR use existing FASTA
        if (params.uniprot_id) {
            // Download and prepare from UniProt using GET_PROTEOME
            GET_PROTEOME(
                output_dir,
                params.reference_store_path,
                params.uniprot_id,
                params.philosopher_reviewed,
                params.philosopher_isoforms,
                params.philosopher_contaminants
            )
            proteome = GET_PROTEOME.out.proteome_fasta
        } else {
            // Use existing FASTA file
            proteome = file(params.reference_proteome)
        }
        
        // NOTE: PHILOSOPHER will be called later after MSFragger produces .pepXML files
        // It combines: annotate, peptideprophet, proteinprophet, filter, report all in one process
        
        // Handle MSFragger config: use provided config file directly, or generate from template
        def msfragger_config
        if (params.msfragger_config_file) {
            // Use provided complete config file directly
            msfragger_config = Channel.value(file(params.msfragger_config_file))
        } else {
            // Generate config from template - only database_name and num_threads will be overridden
            def search_type = params.search_type ?: 'closed'
            def template_file = params.msfragger_config_file ?: 
                               (search_type == "open"         ? "${projectDir}/conf/msfragger/open_fragger.params" :
                                search_type == "nonspecific"  ? "${projectDir}/conf/msfragger/nonspecific_fragger.params" :
                                "${projectDir}/conf/msfragger/closed_fragger.params")
            MSFRAGGER_CONFIG(output_dir, search_type, file(template_file), proteome)
            msfragger_config = MSFRAGGER_CONFIG.out.msfragger_config
        }

        //OPENMS_FILEINFO(output_dir, STAGE_INPUT.out.mzml_files)
        //COMET_SEARCH(output_dir, STAGE_INPUT.out.mzml_files)
        MSFRAGGER(output_dir, STAGE_INPUT.out.mzml_files.map { it[1] }.collect(), msfragger_config, proteome)
        
        // Determine which pepXML files to use for Philosopher based on search type and Crystal-C
        def search_type = params.search_type ?: 'closed'
        
        // Crystal-C mass recalibration (only for open searches)
        if (search_type == "open" && params.run_crystal_c) {
            // Handle Crystal-C config: use provided config file directly, or generate from template
            def crystal_c_config
            if (params.crystal_c_config_file) {
                // Use provided complete config file directly
                crystal_c_config = Channel.value(file(params.crystal_c_config_file))
            } else {
                // Generate config from template
                def template_file = "${projectDir}/conf/crystalc/crystalc.params"
                CRYSTAL_C_CONFIG(output_dir, file(template_file), proteome)
                crystal_c_config = CRYSTAL_C_CONFIG.out.crystal_c_config
            }
            CRYSTAL_C(output_dir, MSFRAGGER.out.pepxml, crystal_c_config)
            philosopher_pepxml = CRYSTAL_C.out.pepxml_recalibrated.collect()
        } else {
            philosopher_pepxml = MSFRAGGER.out.pepxml.collect()
        }
        
        // Run Philosopher (workspace init, annotate, peptideprophet, proteinprophet, filter, report)
        PHILOSOPHER(output_dir, proteome, philosopher_pepxml)
        
        // Run IonQuant for quantification
        IONQUANT(output_dir, PHILOSOPHER.out.psm_tsv, PHILOSOPHER.out.protein_tsv, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())
        
        // Run MSstats for statistical analysis
        MSSTATS(output_dir, IONQUANT.out.msstats, runsheet)

    emit:
        runsheet
        samples
        proteome
}