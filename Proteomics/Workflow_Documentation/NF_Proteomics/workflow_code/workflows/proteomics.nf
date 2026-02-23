include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
//include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { RUNSHEET_TO_MANIFEST } from '../modules/runsheet_to_manifest.nf'
include { STAGE_INPUT } from '../modules/stage_input.nf'
include { RAWBEANS_QC } from '../modules/rawbeans_qc.nf'
include { RAWBEANS_QC_ALL } from '../modules/rawbeans_qc.nf'
include { GET_PROTEOME } from '../modules/get_proteome.nf'
include { CHECK_DECOYS_CONTAMS } from '../modules/check_decoys_contams.nf'
include { MSSTATS } from '../modules/msstats.nf'
include { FRAGPIPE_CONFIG_SETUP } from '../modules/fragpipe_config_setup.nf'
include { FRAGPIPE } from '../modules/fragpipe.nf'
include { RUNSHEET_TO_EXPERIMENT_ANNOTATION } from '../modules/runsheet_to_experiment_annotation.nf'
include { FP_ANALYST } from '../modules/fp_analyst.nf'
include { PMULTIQC } from '../modules/pmultiqc.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


ch_dp_tools_plugin = params.dp_tools_plugin ? 
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file("$projectDir/bin/dp_tools__NF_Proteomics"))

output_dir = Channel.value(file(params.output_dir, type: 'dir', checkIfExists: false))

workflow PROTEOMICS {
    take:
    main:
        // 
        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }
        
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

        // Handle runsheet: generate from ISA or use input params.runsheet file
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
            // Use input params.runsheet file
            runsheet = Channel.value(params.runsheet)
        }

        // Convert the runsheet to a list of tuples of individualsample metadata and input files
        samples = runsheet
            .flatMap { runsheet_path -> 
                samplesheetToList(runsheet_path, "$projectDir/schema_runsheet.json")
            }
            .map { row ->
                def meta = row[0]
                tuple(meta, meta.data_file)
            }

        // Stage input files for each sample
        STAGE_INPUT(output_dir, samples)
        // Run RawBeans QC on each sample's raw data
        RAWBEANS_QC(output_dir, STAGE_INPUT.out.mzml_files)
        // Run RawBeans QC on all samples' raw data
        RAWBEANS_QC_ALL(output_dir, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())

        // Validate database input: cannot specify both uniprot_id and reference_proteome
        if (params.uniprot_id && params.reference_proteome) {
            error "ERROR: Cannot specify both uniprot_id and reference_proteome. Use EITHER uniprot_id to download from UniProt OR reference_proteome to use a custom fasta file."
        }
        
        if (!params.uniprot_id && !params.reference_proteome) {
            error "ERROR: Must specify either uniprot_id or reference_proteome."
        }

        // Download proteome fasta file from UniProt or use input fasta file
        if (params.uniprot_id) {
            // Download and prepare from UniProt using GET_PROTEOME
            GET_PROTEOME(output_dir)
            proteome = GET_PROTEOME.out.proteome_fasta
        } else {
            // Check if fasta file contains decoys and contaminants using input fasta file
            CHECK_DECOYS_CONTAMS(output_dir, file(params.reference_proteome))
            proteome = CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked
        }

        // Generate manifest from runsheet (manifest input disabled for now)
        // if (params.manifest) {
        //     manifest = Channel.fromPath(params.manifest)
        // } else {
        RUNSHEET_TO_MANIFEST(output_dir, runsheet)
        manifest = RUNSHEET_TO_MANIFEST.out.manifest
        // }
        
        ///////////////////////////////////////////////////////////
        // HEADLESS FRAGPIPE PROCESSING:
        ///////////////////////////////////////////////////////////
        
        // https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html
        // https://msfragger-upgrader.nesvilab.org/ionquant/
        // https://msfragger-upgrader.nesvilab.org/diatracer/
        // https://msfragger.arsci.com/upgrader/ (includes ext files)

        // Pass in folder containing externalFragPipe tools
        // |-- tools_folder/
        // |      |-- diaTracer-[version].jar/
        // |      |-- IonQuant-[version].jar/
        // |      |-- MSFragger-[version].jar/
        // |      |-- ext
        // |            |-- bruker/
        // |            |-- thermo/

        ch_fragpipe_tools = params.fragpipe_tools ? Channel.fromPath( params.fragpipe_tools ) : Channel.fromPath("NO_FILE")

        // Point to corresponding FragPipe workflow config file based on params.fragpipe_workflow (string) and optional input params.fragpipe_workflow_config file (path)
        if (params.fragpipe_workflow_config) {
            // Use optional input params.fragpipe_workflow_config file (path)
            fragpipe_config = Channel.fromPath(params.fragpipe_workflow_config)
        } else if (params.fragpipe_workflow) {
            // Or map required input params.fragpipe_workflow (string) to workflow config file path in ${projectDir}/conf/workflows/
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
            error "ERROR: Provide fragpipe_workflow (preset: LFQ-MBR, TMT10, TMT16, TMT16-phospho) and optional fragpipe_workflow_config (path to .workflow file)."
        }

        // Set database file in FragPipe Config file. For TMT workflows, also enable MSstats outputs [STUB]
        FRAGPIPE_CONFIG_SETUP(output_dir, fragpipe_config, proteome)
        FRAGPIPE(output_dir, FRAGPIPE_CONFIG_SETUP.out.fragpipe_config, ch_fragpipe_tools, manifest, proteome, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())

        ///////////////////////////////////////////////////////////
        // END HEADLESS FRAGPIPE
        ///////////////////////////////////////////////////////////

        // Run MultiQC with pmultiqc FragPipe plugin (after FRAGPIPE publishes)
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        ch_fragpipe_output_dir = output_dir
            .combine(FRAGPIPE.out.combined_protein_tsv)
            .map { outdir, _ -> file("${outdir}/FragPipe", type: 'dir') }
        PMULTIQC(output_dir.map { it + "/pmultiqc" }, ch_fragpipe_output_dir)

        // Update experiment_annotation with condition/replicate from runsheet
        RUNSHEET_TO_EXPERIMENT_ANNOTATION(output_dir, FRAGPIPE.out.experiment_annotation, runsheet)

        // Run MSstats (only for LFQ-MBR workflow)
        if (params.fragpipe_workflow == 'LFQ-MBR') {
            // Use the dedicated msstats_csv and fragpipe_manifest outputs from FRAGPIPE
            MSSTATS(output_dir, runsheet, FRAGPIPE.out.fragpipe_manifest, FRAGPIPE.out.msstats_csv)

            // Run FragPipe-Analyst R on combined_protein.tsv and combined_peptide.tsv
            ch_fp_analyst_protein = FRAGPIPE.out.combined_protein_tsv
                .combine(RUNSHEET_TO_EXPERIMENT_ANNOTATION.out.experiment_annotation)
                .map { quant_file, exp_anno -> tuple("protein", quant_file, exp_anno) }
            ch_fp_analyst_peptide = FRAGPIPE.out.combined_peptide_tsv
                .combine(RUNSHEET_TO_EXPERIMENT_ANNOTATION.out.experiment_annotation)
                .map { quant_file, exp_anno -> tuple("peptide", quant_file, exp_anno) }
            ch_fp_analyst_inputs = ch_fp_analyst_protein.mix(ch_fp_analyst_peptide)

            FP_ANALYST(output_dir, ch_fp_analyst_inputs)
        }

        // Software version capturing
        ch_software_versions = Channel.empty()
            | mix(STAGE_INPUT.out.versions)
            | mix(RAWBEANS_QC_ALL.out.versions)
            | mix(FRAGPIPE.out.versions)
            | mix(PMULTIQC.out.versions)
            | mix(RUNSHEET_TO_MANIFEST.out.versions)
            | mix(MSSTATS.out.versions)
            | mix(FP_ANALYST.out.versions)
        ch_software_versions
            | unique
            | collectFile(newLine: true)
            | set { ch_final_software_versions }
        SOFTWARE_VERSIONS(output_dir, ch_final_software_versions)

    emit:
        runsheet
        samples
        proteome
}