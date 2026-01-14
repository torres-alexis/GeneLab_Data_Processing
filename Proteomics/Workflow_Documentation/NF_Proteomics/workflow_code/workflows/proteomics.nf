include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
//include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { RUNSHEET_TO_MANIFEST } from '../modules/runsheet_to_manifest.nf'
include { STAGE_INPUT } from '../modules/stage_input.nf'
include { RAWBEANS_QC } from '../modules/rawbeans_qc.nf'
include { RAWBEANS_QC_ALL } from '../modules/rawbeans_qc.nf'
//include { OPENMS_FILEINFO } from '../modules/openms_fileinfo.nf'
//include { MSFRAGGER } from '../modules/msfragger.nf'
//include { MSFRAGGER_CONFIG } from '../modules/msfragger.nf'
include { GET_PROTEOME } from '../modules/get_proteome.nf'
//include { CRYSTAL_C } from '../modules/crystal_c.nf'
//include { CRYSTAL_C_CONFIG } from '../modules/crystal_c.nf'
//include { PHILOSOPHER } from '../modules/philosopher.nf'
//include { IONQUANT } from '../modules/ionquant.nf'
//include { MSSTATS } from '../modules/msstats.nf'
include { FRAGPIPE_CONFIG_SETUP } from '../modules/fragpipe_config_setup.nf'
include { FRAGPIPE } from '../modules/fragpipe.nf'
//include { PEPXML_TO_MZID } from '../modules/pepxml_to_mzid.nf'
//include { PMULTIQC } from '../modules/pmultiqc.nf'



include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

import groovy.json.JsonSlurper

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
        if ( params.runsheet == null ) {
            // Generate runsheet from ISA archive // STUB //
            if ( params.isa_archive == null ) {
                FETCH_ISA( output_dir, osd_accession, glds_accession )
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, FETCH_ISA.out.isa_archive, ch_dp_tools_plugin ) // STUB //
                runsheet = ISA_TO_RUNSHEET.out.runsheet.map { it.toString() }
            } else {
                ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, params.isa_archive, ch_dp_tools_plugin ) // STUB //
                runsheet = ISA_TO_RUNSHEET.out.runsheet.map { it.toString() }
            }
        } else {
            // Use provided runsheet
            runsheet = Channel.value(params.runsheet)
        }

        // Convert the runsheet to a list of samples and run processing
        samples = runsheet
            .flatMap { runsheet_path -> 
                samplesheetToList(runsheet_path, "$projectDir/schema_runsheet.json")
            }
            .map { meta, input_file -> 
                tuple(meta, input_file)
            }
        
        // samples.view { meta, input_file ->
        //     """Sample: ${meta.id}, Input File: ${input_file}"""
        // }

        // Stage input files for each sample
        STAGE_INPUT(output_dir, samples)
        // Run QC on each sample's raw data
        RAWBEANS_QC(output_dir, STAGE_INPUT.out.mzml_files)
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

        // Generate manifest if not supplied
        if (params.manifest) {
            manifest = Channel.fromPath(params.manifest)
        } else {
            RUNSHEET_TO_MANIFEST(output_dir, runsheet)
            manifest = RUNSHEET_TO_MANIFEST.out.manifest
        }
        
        ///////////////////////////////////////////////////////////
        // HEADLESS FRAGPIPE PROCESSING:
        ///////////////////////////////////////////////////////////
        ch_fragpipe_tools = params.fragpipe_tools ? Channel.fromPath( params.fragpipe_tools ) : Channel.fromPath("NO_FILE")

        // Determine fragpipe_config: use fragpipe_workflow_config if provided (takes precedence), otherwise map fragpipe_workflow string to workflow file
        if (params.fragpipe_workflow_config) {
            // Use provided config file (takes precedence over workflow string)
            fragpipe_config = Channel.fromPath(params.fragpipe_workflow_config)
        } else if (params.fragpipe_workflow) {
            // Map fragpipe_workflow string to workflow file path
            def workflow_file
            switch(params.fragpipe_workflow) {
                case 'TMT10':
                    workflow_file = "${projectDir}/conf/workflows/TMT10.workflow"
                    break
                case 'TMT16':
                    workflow_file = "${projectDir}/conf/workflows/TMT16.workflow"
                    break
                case 'TMT16-phospho':
                    workflow_file = "${projectDir}/conf/workflows/TMT16-phospho.workflow"
                    break
                case 'LFQ-MBR':
                default:
                    workflow_file = "${projectDir}/conf/workflows/LFQ-MBR.workflow"
                    break
            }
            fragpipe_config = Channel.value(file(workflow_file))
        } else {
            error "ERROR: Either fragpipe_workflow_config or fragpipe_workflow parameter must be provided. Use fragpipe_workflow_config to specify a custom workflow file, or fragpipe_workflow to use a preset (LFQ-MBR, TMT10, TMT16, TMT16-phospho)."
        }

        FRAGPIPE_CONFIG_SETUP(output_dir, fragpipe_config, proteome)
        
        slurp = new JsonSlurper()

        // View enabled workflow modules
        FRAGPIPE_CONFIG_SETUP.out.fragpipe_json
            .map { json_str ->
                def config = slurp.parseText(json_str)
                
                // Custom display order - matches FragPipe execution order from log
                // Workflow order columns: [LFQ-MBR | TMT-10 | TMT-16 | TMT-16-phospho]
                // Execution order: MSFragger -> MSBooster -> Validation -> Percolator -> ProteinProphet -> IonQuant -> Reporting
                def display_order = [
                    'msfragger.run-msfragger',           // 1 | 1 | 1 | 1 - Search engine
                    'msbooster.run-msbooster',           // 2 | 2 | 2 | 2 - MS2 prediction/boosting
                    'run-validation-tab',                 // Flag - Validation tab
                    'run-psm-validation',                // Flag - Enables PSM validation
                    'percolator.run-percolator',         // 3 | 3 | 3 | 3 - PSM validation tool
                    'protein-prophet.run-protein-prophet', // 4 | 4 | 4 | 4 - Protein inference
                    'peptide-prophet.run-peptide-prophet', // (if enabled, runs after Percolator)
                    'ptmprophet.run-ptmprophet',         // (if enabled, runs for PTM localization)
                    'phi-report.run-report',             // 5 | 5 | 5 | 5 - Generate reports (runs after ProteinProphet, before IonQuant)
                    'quantitation.run-label-free-quant', // Flag - Configures IonQuant for label-free quant
                    'ionquant.run-ionquant',             // 6 | 6 | 6 | 6 - Quantification tool
                    // Other modules (not in LFQ-MBR workflow)
                    'diann.run-dia-nn',
                    'diaumpire.run-diaumpire',
                    'diatracer.run-diatracer',
                    'freequant.run-freequant',
                    'tmtintegrator.run-tmtintegrator',
                    'crystalc.run-crystalc',
                    'fpop.fragpipe.fpop.run-fpop',
                    'mbg.run-mbg',
                    'metaproteomics.run-metaproteomics',
                    'opair.run-opair',
                    'saintexpress.run-saint-express',
                    'skyline.run-skyline',
                    'speclibgen.run-speclibgen',
                    'transfer-learning.run-transfer-learning'
                ]
                
                // Filter for enabled modules (value == "true") in display order
                def enabled_modules = display_order.findAll { flag ->
                    config[flag] == "true"
                }
                
                // Format output
                if (enabled_modules.size() > 0) {
                    return "Enabled FragPipe modules:\n" + enabled_modules.collect { "  ✓ ${it}" }.join("\n")
                } else {
                    return "No modules enabled"
                }
            }
            .view { it }


        FRAGPIPE(output_dir, FRAGPIPE_CONFIG_SETUP.out.fragpipe_config, ch_fragpipe_tools, manifest, proteome, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())

        // // Run MultiQC with pmultiqc fragpipe plugin
        // ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        // ch_fragpipe_output_dir = output_dir.map { file("${it}/FragPipe", type: 'dir') }
        // PMULTIQC( output_dir.map { it + "/pmultiqc" }, ch_fragpipe_output_dir)

        ///////////////////////////////////////////////////////////
        // END HEADLESS FRAGPIPE
        ///////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////
        // MANUAL FRAGPIPE PROCESSING:
        ///////////////////////////////////////////////////////////

        // // Handle MSFragger config: use provided config file directly, or generate from template
        // def msfragger_config
        // if (params.msfragger_config_file) {
        //     // Use provided complete config file directly
        //     msfragger_config = Channel.value(file(params.msfragger_config_file))
        // } else {
        //     // Generate config from template - only database_name and num_threads will be overridden
        //     def search_type = params.search_type ?: 'closed'
        //     def template_file = params.msfragger_config_file ?: 
        //                        (search_type == "open"         ? "${projectDir}/conf/msfragger/open_fragger.params" :
        //                         search_type == "nonspecific"  ? "${projectDir}/conf/msfragger/nonspecific_fragger.params" :
        //                         "${projectDir}/conf/msfragger/closed_fragger.params")
        //     MSFRAGGER_CONFIG(output_dir, search_type, file(template_file), proteome)
        //     msfragger_config = MSFRAGGER_CONFIG.out.msfragger_config
        // }

        // //OPENMS_FILEINFO(output_dir, STAGE_INPUT.out.mzml_files)
        // MSFRAGGER(output_dir, STAGE_INPUT.out.mzml_files.map { it[1] }.collect(), msfragger_config, proteome)
        
        // // Determine which pepXML files to use for Philosopher based on search type and Crystal-C
        // def search_type = params.search_type ?: 'closed'
        
        // // Crystal-C mass recalibration (only for open searches)
        // if (search_type == "open" && params.run_crystal_c) {
        //     // Handle Crystal-C config: use provided config file directly, or generate from template
        //     def crystal_c_config
        //     if (params.crystal_c_config_file) {
        //         // Use provided complete config file directly
        //         crystal_c_config = Channel.value(file(params.crystal_c_config_file))
        //     } else {
        //         // Generate config from template
        //         def template_file = "${projectDir}/conf/crystalc/crystalc.params"
        //         CRYSTAL_C_CONFIG(output_dir, file(template_file), proteome)
        //         crystal_c_config = CRYSTAL_C_CONFIG.out.crystal_c_config
        //     }
        //     CRYSTAL_C(output_dir, MSFRAGGER.out.pepxml, crystal_c_config)
        //     philosopher_pepxml = CRYSTAL_C.out.pepxml_recalibrated.collect()
        // } else {
        //     philosopher_pepxml = MSFRAGGER.out.pepxml.collect()
        // }
        
        // // Run Philosopher (workspace init, annotate, peptideprophet, proteinprophet, filter, report)
        // PHILOSOPHER(output_dir, proteome, philosopher_pepxml)
        
        // // Run IonQuant for quantification
        // IONQUANT(output_dir, PHILOSOPHER.out.psm_tsv, PHILOSOPHER.out.protein_tsv, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())
        
        ///////////////////////////////////////////////////////////
        // END MANUAL FRAGPIPE PROCESSING
        ///////////////////////////////////////////////////////////

        // Run MSstats for statistical analysis (only for LFQ-MBR workflow)
        if (params.fragpipe_workflow == 'LFQ-MBR') {
            // Use the dedicated msstats_csv and fragpipe_manifest outputs from FRAGPIPE
            MSSTATS(output_dir, runsheet, FRAGPIPE.out.fragpipe_manifest, FRAGPIPE.out.msstats_csv)
        }


    emit:
        runsheet
        samples
        proteome
}