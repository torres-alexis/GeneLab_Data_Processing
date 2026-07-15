include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { GET_PROTEOME } from '../modules/get_proteome.nf'
include { FETCH_REFERENCE_PROTEOME } from '../modules/fetch_reference_proteome.nf'
include { CHECK_DECOYS_CONTAMS } from '../modules/check_decoys_contams.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { ISA_TO_TMT_SHEETS } from '../modules/isa_to_tmt_sheets.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { STAGE_INPUT } from '../modules/stage_input.nf'
// include { RAWBEANS_QC } from '../modules/rawbeans_qc.nf'
include { RAWBEANS_QC_ALL } from '../modules/rawbeans_qc.nf'
include { FRAGPIPE_CONFIG_SETUP } from '../modules/fragpipe_config_setup.nf'
include { FRAGPIPE_METADATA_SETUP } from '../modules/fragpipe_metadata_setup.nf'
include { FRAGPIPE } from '../modules/fragpipe.nf'
include { CLEAN_FRAGPIPE_TABLES } from '../modules/clean_fragpipe_tables.nf'
include { ZIP_FRAGPIPE_OUTPUTS } from '../modules/zip_fragpipe_outputs.nf'
include { PMULTIQC } from '../modules/pmultiqc.nf'
include { MSSTATS } from '../modules/msstats.nf'
include { MSSTATSTMT } from '../modules/msstatstmt.nf'
include { FRAGPIPEANALYSTR }  from '../modules/fragpipeanalystr.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROCESSED_PROTOCOL } from '../modules/generate_protocol.nf'

include { FILTER_TECH_REPS } from '../modules/filter_tech_reps.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow PROTEOMICS {
    main:
        ch_dp_tools_plugin = params.dp_tools_plugin ?
            Channel.value(file(params.dp_tools_plugin)) :
            Channel.value(file("$projectDir/bin/dp_tools__NF_Proteomics"))

        output_dir = Channel.value(file(params.output_dir, type: 'dir', checkIfExists: false))

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

        // Emission from map/combine is queue → use .first() to get first value
        ch_out_dir = output_dir.first()

        // TMT: data_sheet + sample_sheet (runsheet not used). LFQ: runsheet (or generate from ISA).
        def is_tmt = params.fragpipe_workflow?.startsWith('TMT')
        def sheet

        if ( is_tmt ) {
            if ( params.data_sheet && params.sample_sheet ) {
                sheet = Channel.fromPath(params.data_sheet, checkIfExists: true)
            // } else if ( params.accession ) {
            //     // STUB: ISA to TMT sheets. Implement dpt-isa-to-tmt-sheets in dp_tools later.
            //     if ( params.isa_archive == null ) {
            //         FETCH_ISA( output_dir, osd_accession, glds_accession )
            //         ISA_TO_TMT_SHEETS( output_dir, osd_accession, glds_accession, FETCH_ISA.out.isa_archive, ch_dp_tools_plugin )
            //     } else {
            //         ISA_TO_TMT_SHEETS( output_dir, osd_accession, glds_accession, file(params.isa_archive), ch_dp_tools_plugin )
            //     }
            //     sheet = ISA_TO_TMT_SHEETS.out.data_sheet.map { it.toString() }
            } else {
                error "TMT workflows require --data_sheet and --sample_sheet."
            }
        } else {
            // LFQ: runsheet from params or ISA
            if ( params.runsheet == null ) {
                if ( !params.accession ) {
                    error "When generating runsheet from ISA, --accession (OSD-### or GLDS-###) is required."
                }
                if ( params.isa_archive == null ) {
                    FETCH_ISA( output_dir, osd_accession, glds_accession )
                    ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, FETCH_ISA.out.isa_archive, ch_dp_tools_plugin )
                    sheet = ISA_TO_RUNSHEET.out.runsheet.map { it.toString() }
                } else {
                    ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, file(params.isa_archive), ch_dp_tools_plugin )
                    sheet = ISA_TO_RUNSHEET.out.runsheet.map { it.toString() }
                }
            } else {
                sheet = Channel.value(params.runsheet)
            }
        }

        def sheet_in = sheet
        def sheet_fp
        if (params.first_technical_replicate_only) {
            FILTER_TECH_REPS(ch_out_dir, sheet_in)
            sheet_fp = FILTER_TECH_REPS.out.filtered
        } else {
            sheet_fp = sheet_in
        }

        // Convert sheet to list of samples (file metadata + path)
        def sheet_schema = is_tmt ? "$projectDir/schema_data_sheet.json" : "$projectDir/schema_runsheet.json"
        // Stage + RawBeans QC on full runsheet (all technical replicates)
        samples_full = sheet_in
            .flatMap { sheet_path ->
                samplesheetToList(sheet_path, sheet_schema)
            }
            .map { row -> row[0] }
        // FragPipe manifest + search use first tech rep only when params.first_technical_replicate_only
        def samples_fp
        if (params.first_technical_replicate_only) {
            samples_fp = sheet_fp
                .flatMap { sheet_path ->
                    samplesheetToList(sheet_path, sheet_schema)
                }
                .map { row -> row[0] }
        } else {
            samples_fp = samples_full
        }

        // TMT: parse sample_sheet for organism etc. (sample-centric; same meta structure as runsheet)
        sample_sheet_rows = (is_tmt && params.sample_sheet)
            ? Channel.fromPath(params.sample_sheet).flatMap { samplesheetToList(it, "$projectDir/schema_sample_sheet.json") }
            : Channel.empty()

        // Stage input mzML files for each sample / fraction
        STAGE_INPUT(ch_out_dir, samples_full)
        // Run RawBeans QC on each sample's raw data
        // RAWBEANS_QC(ch_out_dir, STAGE_INPUT.out.mzml_files)
        // Run RawBeans QC on all samples' raw data
        RAWBEANS_QC_ALL(output_dir, STAGE_INPUT.out.mzml_files.map { it[1] }.collect())

        // Organism from runsheet / sample_sheet (one per dataset). Used for reference table lookups.
        ch_organism_sci = (is_tmt ? sample_sheet_rows : samples_fp)
            | first
            | map { row ->
                def meta = row instanceof Map ? row : row[0]
                meta?.organism ? meta.organism.toString().replaceAll(" ","_").toLowerCase().trim() : ""
            }

        if (params.reference_table) {
            PARSE_ANNOTATIONS_TABLE(Channel.value(params.reference_table), ch_organism_sci)
        }

        // Reference proteome: explicit params override; else resolve from GL-DPPD-7110-A table + runsheet organism.
        if (params.uniprot_id && params.reference_proteome) {
            error "ERROR: Cannot specify both uniprot_id and reference_proteome."
        }

        if (params.uniprot_id) {
            GET_PROTEOME(output_dir)
            proteome = GET_PROTEOME.out.proteome_fasta
        } else if (params.reference_proteome) {
            CHECK_DECOYS_CONTAMS(output_dir, file(params.reference_proteome))
            proteome = CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked
        } else if (params.reference_table) {
            ch_proteome_source = PARSE_ANNOTATIONS_TABLE.out.proteome_source
                .map { src ->
                    if (!src) {
                        error "ERROR: No proteome in reference_table for this organism. Specify uniprot_id or reference_proteome."
                    }
                    src
                }
            FETCH_REFERENCE_PROTEOME(output_dir, ch_proteome_source)
            CHECK_DECOYS_CONTAMS(output_dir, FETCH_REFERENCE_PROTEOME.out.proteome_fasta)
            proteome = CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked
        } else {
            error "ERROR: Must specify uniprot_id, reference_proteome, or reference_table with a matching organism row."
        }

        // Generate manifest: TMT [data_sheet, sample_sheet]; LFQ [runsheet]
        ch_sample_sheet = (is_tmt && params.sample_sheet) ? Channel.fromPath(params.sample_sheet) : Channel.value(file("${projectDir}/bin/placeholder"))
        ch_sheets = sheet_fp.combine(ch_sample_sheet).map { s, ss -> [s, ss] }
        FRAGPIPE_METADATA_SETUP(output_dir, ch_sheets)
        manifest = FRAGPIPE_METADATA_SETUP.out.manifest
        // }
        
        ///////////////////////////////////////////////////////////
        // HEADLESS FRAGPIPE PROCESSING:
        ///////////////////////////////////////////////////////////
        
        // https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html
        // https://msfragger-upgrader.nesvilab.org/ionquant/
        // https://msfragger-upgrader.nesvilab.org/diatracer/
        // https://msfragger.arsci.com/upgrader/ (includes ext files)

        // Pass in folder containing external FragPipe tools
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
            // v2 parser: classic switch/break is not accepted here; keep if-chain for preset mapping
            if (params.fragpipe_workflow == 'TMT10') {
                workflow_file = "${projectDir}/conf/workflows/TMT10.workflow"
            } else if (params.fragpipe_workflow == 'TMT16') {
                workflow_file = "${projectDir}/conf/workflows/TMT16.workflow"
            } else if (params.fragpipe_workflow == 'TMT16-phospho') {
                workflow_file = "${projectDir}/conf/workflows/TMT16-phospho.workflow"
            } else {
                workflow_file = "${projectDir}/conf/workflows/LFQ-MBR.workflow"
            }
            fragpipe_config = Channel.value(file(workflow_file))
        } else {
            error "ERROR: Provide fragpipe_workflow (preset: LFQ-MBR, TMT10, TMT16, TMT16-phospho) and optional fragpipe_workflow_config (path to .workflow file)."
        }

        // Set database file in FragPipe Config file. For TMT workflows, also enable MSstats outputs [STUB]
        FRAGPIPE_CONFIG_SETUP(output_dir, fragpipe_config, proteome)
        // Join on id — do not combine(collect()) id lists: Groovy flattens List into tuple slots
        ch_fragpipe_mzml = params.first_technical_replicate_only
            ? STAGE_INPUT.out.mzml_files
                .map { meta, mzml -> tuple(meta.id.toString(), mzml) }
                .join(samples_fp.map { m -> tuple(m.id.toString(), 1) }, by: 0)
                .map { id, mzml, fpMarker -> mzml }
                .collect()
            : STAGE_INPUT.out.mzml_files.map { it[1] }.collect()

        FRAGPIPE(output_dir, FRAGPIPE_CONFIG_SETUP.out.fragpipe_config, ch_fragpipe_tools, manifest, proteome, ch_fragpipe_mzml, FRAGPIPE_METADATA_SETUP.out.experiment_annotation)

        ///////////////////////////////////////////////////////////
        // END HEADLESS FRAGPIPE
        ///////////////////////////////////////////////////////////

        // Run pmultiqc with FragPipe plugin
        ch_fragpipe_output_dir = FRAGPIPE.out.fragpipe_manifest.map { fragpipe_manifest -> fragpipe_manifest.parent }
        PMULTIQC(output_dir.map { it + "/pmultiqc" }, ch_fragpipe_output_dir)

        // Pass in FragPipe tables to clean and publish; run downstream modules with original tables
        ch_fragpipe_tables_to_clean = (params.fragpipe_workflow?.startsWith('TMT') ?
            FRAGPIPE.out.tmt_report_tables :
            FRAGPIPE.out.combined_tables
        ).collect()
        CLEAN_FRAGPIPE_TABLES(
            ch_out_dir,
            ch_fragpipe_tables_to_clean,
            FRAGPIPE_METADATA_SETUP.out.experiment_annotation
        )

        ZIP_FRAGPIPE_OUTPUTS(ch_out_dir, ch_fragpipe_output_dir)

        // Resolve gene_annotations_url for DE_results: direct file, or organism + reference table
        gene_annotations_url = Channel.value(null)
        if (params.gene_annotations_file) {
            gene_annotations_url = Channel.value(params.gene_annotations_file)
        } else if (params.reference_table) {
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }
        ch_lfq_versions = Channel.empty()
        ch_tmt_versions = Channel.empty()
        if (params.fragpipe_workflow == 'LFQ-MBR') {
            MSSTATS(output_dir, FRAGPIPE_METADATA_SETUP.out.experiment_annotation, FRAGPIPE.out.msstats_csv)
            ch_lfq_versions = MSSTATS.out.versions
        }
        if (params.fragpipe_workflow?.startsWith('TMT')) {
            MSSTATSTMT(
                output_dir,
                FRAGPIPE_METADATA_SETUP.out.msstats_tmt_annotation,
                FRAGPIPE.out.msstats_csv
            )
            ch_tmt_versions = MSSTATSTMT.out.versions
        }

        // FRAGPIPEANALYSTR (FragPipeAnalystR): levels from params or workflow default. TMT uses tmt-report abundance/ratio; LFQ uses combined_*.
        def fp_levels = params.fp_analyst_levels ?
            params.fp_analyst_levels.toString().split(',')*.trim().findAll { it } :
            (params.fragpipe_workflow?.startsWith('TMT') ? ['protein','gene','peptide','site'] : ['protein','peptide'])
        def tmt_quant = (params.fp_analyst_tmt_quant_type ?: 'abundance').toLowerCase().startsWith('ratio') ? 'ratio' : 'abundance'
        def takeFirst = { p -> p instanceof List ? p[0] : p }

        ch_fp_analyst_inputs = Channel.empty()
        if (fp_levels.contains('protein')) {
            def ch_protein = params.fragpipe_workflow?.startsWith('TMT') ?
                (tmt_quant == 'ratio' ? FRAGPIPE.out.ratio_protein : FRAGPIPE.out.abundance_protein).map(takeFirst) :
                FRAGPIPE.out.combined_protein
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_protein.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("protein", q, e) })
        }
        if (fp_levels.contains('peptide')) {
            def ch_peptide = params.fragpipe_workflow?.startsWith('TMT') ?
                (tmt_quant == 'ratio' ? FRAGPIPE.out.ratio_peptide : FRAGPIPE.out.abundance_peptide).map(takeFirst) :
                FRAGPIPE.out.combined_peptide
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_peptide.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("peptide", q, e) })
        }
        if (params.fragpipe_workflow?.startsWith('TMT') && fp_levels.contains('gene')) {
            def ch_gene = (tmt_quant == 'ratio' ? FRAGPIPE.out.ratio_gene : FRAGPIPE.out.abundance_gene).map(takeFirst)
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_gene.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("gene", q, e) })
        }
        if (params.fragpipe_workflow?.startsWith('TMT') && fp_levels.contains('site')) {
            def ch_site = (tmt_quant == 'ratio' ? FRAGPIPE.out.ratio_single_site : FRAGPIPE.out.abundance_single_site).map(takeFirst)
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_site.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("site", q, e) })
        }
        ch_fp_with_annot = ch_fp_analyst_inputs.combine(gene_annotations_url)
        FRAGPIPEANALYSTR(ch_out_dir, ch_fp_with_annot)

        ch_lfq_versions = ch_lfq_versions.mix(ch_tmt_versions).mix(FRAGPIPEANALYSTR.out.versions)

        // Software version capturing
        ch_software_versions = Channel.empty()
            | mix(STAGE_INPUT.out.versions)
            | mix(RAWBEANS_QC_ALL.out.versions)
            | mix(FRAGPIPE.out.versions)
            | mix(PMULTIQC.out.versions)
            | mix(FRAGPIPE_METADATA_SETUP.out.versions)
            | mix(ch_lfq_versions)
            | mix(ch_tmt_versions)
        ch_software_versions
            | unique
            | collectFile(newLine: true)
            | set { ch_final_software_versions }
        SOFTWARE_VERSIONS(output_dir, ch_final_software_versions)

        if (params.uniprot_id) {
            ch_protocol_uniprot = Channel.value(params.uniprot_id.toString())
        } else if (params.reference_table) {
            ch_protocol_uniprot = PARSE_ANNOTATIONS_TABLE.out.uniprot_id.map { uid -> uid?.toString() ?: '' }
        } else {
            ch_protocol_uniprot = Channel.value('')
        }
        ch_protocol_reference_table = Channel.value(params.reference_table?.toString() ?: '')
        ch_protocol_reference_proteome = Channel.value(params.reference_proteome?.toString() ?: '')
        GENERATE_PROCESSED_PROTOCOL(
            ch_out_dir,
            SOFTWARE_VERSIONS.out.software_versions,
            proteome,
            ch_protocol_uniprot,
            ch_protocol_reference_table,
            ch_protocol_reference_proteome
        )
}