include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { GET_PROTEOME } from '../modules/get_proteome.nf'
include { FETCH_REFERENCE_PROTEOME } from '../modules/fetch_reference_proteome.nf'
include { CHECK_DECOYS_CONTAMS } from '../modules/check_decoys_contams.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { COPY_INPUT } from '../modules/copy_input.nf'
include { FETCH_INPUT } from '../modules/fetch_input.nf'
include { RAWBEANS_QC_ALL } from '../modules/rawbeans_qc.nf'
include { FRAGPIPE_CONFIG_SETUP } from '../modules/fragpipe_config_setup.nf'
include { FRAGPIPE_METADATA_SETUP } from '../modules/fragpipe_metadata_setup.nf'
include { FRAGPIPE } from '../modules/fragpipe.nf'
include { CLEAN_FRAGPIPE_TABLES } from '../modules/clean_fragpipe_tables.nf'
include { ZIP_FRAGPIPE_OUTPUTS } from '../modules/zip_fragpipe_outputs.nf'
include { PMULTIQC } from '../modules/pmultiqc.nf'
include { DROP_DECOYS_CONTAMS } from '../modules/drop_decoys_contams.nf'
include { MSSTATS } from '../modules/msstats.nf'
include { MSSTATSTMT } from '../modules/msstatstmt.nf'
include { FRAGPIPEANALYSTR }  from '../modules/fragpipeanalystr.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROCESSED_PROTOCOL } from '../modules/generate_protocol.nf'
include { FILTER_TECH_REPS } from '../modules/filter_tech_reps.nf'
include { VV_STEP; VV_CONCAT_FILTER } from '../modules/vv_step.nf'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

def as_list(x) {
    if (x == null) return []
    if (x instanceof Collection && !(x instanceof CharSequence)) return x as List
    return [x]
}

def files_only(items) {
    return items.collectMany { x -> x instanceof Map ? [] : as_list(x) }
}

def dest_file(root, dir, f) {
    return ["${root}/${dir}/${f.name}".toString(), f]
}

def rel_under(f, marker) {
    def s = f.toString().replace('\\', '/')
    def key = "/${marker}/"
    def i = s.lastIndexOf(key)
    return i >= 0 ? s.substring(i + key.length()) : f.name
}

def pub(ch, root, dir) {
    return ch.combine(root).combine(channel.value(dir)).flatMap { row ->
        def items = as_list(row)
        def d = items[-1]
        def r = items[-2]
        files_only(items[0..-3]).collect { f -> dest_file(r, d, f) }
    }
}

def pub_fpar(ch, root) {
    return ch.combine(root).flatMap { data_type, files, r ->
        as_list(files).collect { f ->
            def rel = rel_under(f, 'output')
            def dest = rel.toLowerCase().endsWith('.rdata')
                ? "${r}/processing_info/${rel}".toString()
                : "${r}/FragPipeAnalystR/${data_type}/${rel}".toString()
            [dest, f]
        }
    }
}

def takeFirst(p) {
    return p instanceof List ? p[0] : p
}

def isRemoteSrc(meta) {
    def src = meta.data_file.toString()
    return src.startsWith('s3://') || src.startsWith('http://') || src.startsWith('https://')
}

def remapLfqFractionIds(metas, is_tmt) {
    if (is_tmt) {
        return metas
    }
    def counts = metas.countBy { it.id?.toString() }
    return metas.collect { meta ->
        def sid = meta.id?.toString()
        if (!sid || (counts[sid] ?: 0) <= 1) {
            return meta
        }
        def base = file(meta.data_file.toString()).getName()
        def stem = base.replaceAll(/(?i)\.mzML$/, '')
        return meta + [id: stem]
    }
}

def existing_under(root, patterns) {
    return patterns.collectMany { pat ->
        def p = "${root}/${pat}".toString()
        def hits
        if (p.contains('*') || p.contains('?')) {
            hits = files(p, checkIfExists: false)
        } else {
            hits = [file(p)]
        }
        as_list(hits).findAll { f -> f.exists() && f.isFile() }
    }.unique { it.toString() }
}

def require_under(root, patterns, label) {
    def hits = existing_under(root, patterns)
    if (!hits) {
        error "entry_point=fragpipe_output: no ${label} under ${root} (tried ${patterns.join(', ')})"
    }
    return hits
}

def resolve_fp_root(raw) {
    def d = file(raw)
    if (!d.exists()) {
        error "entry_point=fragpipe_output: ${raw} does not exist."
    }
    if (!d.isDirectory()) {
        error "entry_point=fragpipe_output: ${raw} is not a directory."
    }
    def candidates = [d, file("${d}/output"), file("${d}/FragPipe")]
    def hit = candidates.find { c -> c.exists() && c.isDirectory() && file("${c}/msstats.csv").exists() }
    if (!hit) {
        error "entry_point=fragpipe_output: no msstats.csv in ${raw}."
    }
    return hit
}

def is_fp_workdir(root) {
    def kids = root.listFiles()
    return file("${root}/psm.tsv").exists() || (kids != null && kids.any { f -> f.name.startsWith('log_') && f.name.endsWith('.txt') })
}

def files_ch(files) {
    if (!files) {
        return Channel.empty()
    }
    return Channel.fromList(files)
}

def tmt_table_pats(kind, quant) {
    def prefix = quant == 'ratio' ? 'ratio' : 'abundance'
    return ["tmt-report/${prefix}_${kind}_*.tsv", "${prefix}_${kind}_*.tsv"]
}

workflow PROTEOMICS {
    main:
        def ep = (params.entry_point ?: 'mzml').toString().trim().toLowerCase()
        if (!(ep in ['mzml', 'fragpipe_output'])) {
            error "params.entry_point must be mzml or fragpipe_output."
        }
        if (ep == 'fragpipe_output') {
            if (!params.fragpipe_output) {
                error "--fragpipe_output is required when entry_point=fragpipe_output."
            }
        } else if (params.fragpipe_output) {
            error "--fragpipe_output is only used with --entry_point fragpipe_output."
        }
        println "Entry point: '${ep}'"

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
        if (params.results_dir) {
            ch_root = Channel.value(params.results_dir.toString())
        } else if (params.accession) {
            ch_root = osd_accession
        } else {
            ch_root = Channel.value('results')
        }
        ch_published = Channel.empty()

        // TMT: data_sheet + sample_sheet (runsheet not used). LFQ: runsheet (or generate from ISA).
        def is_tmt = params.fragpipe_workflow?.startsWith('TMT')
        def sheet

        if ( is_tmt ) {
            if ( params.data_sheet && params.sample_sheet ) {
                sheet = Channel.fromPath(params.data_sheet, checkIfExists: true)
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
                    ch_published = ch_published
                        .mix(pub(FETCH_ISA.out.isa_archive, ch_root, 'Metadata'))
                        .mix(pub(ISA_TO_RUNSHEET.out.runsheet, ch_root, 'Metadata'))
                        .mix(pub(ISA_TO_RUNSHEET.out.isa_copy, ch_root, 'Metadata'))
                } else {
                    ISA_TO_RUNSHEET( output_dir, osd_accession, glds_accession, file(params.isa_archive), ch_dp_tools_plugin )
                    sheet = ISA_TO_RUNSHEET.out.runsheet.map { it.toString() }
                    ch_published = ch_published
                        .mix(pub(ISA_TO_RUNSHEET.out.runsheet, ch_root, 'Metadata'))
                        .mix(pub(ISA_TO_RUNSHEET.out.isa_copy, ch_root, 'Metadata'))
                }
            } else {
                sheet = Channel.value(params.runsheet)
            }
        }

        def tech_rep = (params.tech_rep ?: 'first').toString().trim().toLowerCase()
        if (!(tech_rep in ['first', 'all'])) {
            error "params.tech_rep must be first or all."
        }
        def keep_first_tr = (tech_rep == 'first')

        def sheet_in = sheet
        def sheet_fp
        if (keep_first_tr) {
            FILTER_TECH_REPS(ch_out_dir, sheet_in)
            sheet_fp = FILTER_TECH_REPS.out.filtered
            ch_published = ch_published
                .mix(pub(FILTER_TECH_REPS.out.original, ch_root, 'Metadata'))
                .mix(pub(FILTER_TECH_REPS.out.filtered_publish, ch_root, 'Metadata'))
                .mix(pub(FILTER_TECH_REPS.out.drop_log, ch_root, 'Metadata'))
        } else {
            sheet_fp = sheet_in
        }

        def sheet_schema = is_tmt ? "$projectDir/schema_data_sheet.json" : "$projectDir/schema_runsheet.json"
        samples_full = sheet_in
            .flatMap { sheet_path ->
                samplesheetToList(sheet_path, sheet_schema)
            }
            .map { row -> row[0] }
            .toList()
            .flatMap { metas -> remapLfqFractionIds(metas, is_tmt) }
        // FragPipe manifest + search use collapsed sheet when tech_rep=first
        def samples_fp
        if (keep_first_tr) {
            samples_fp = sheet_fp
                .flatMap { sheet_path ->
                    samplesheetToList(sheet_path, sheet_schema)
                }
                .map { row -> row[0] }
                .toList()
                .flatMap { metas -> remapLfqFractionIds(metas, is_tmt) }
        } else {
            samples_fp = samples_full
        }

        // TMT: parse sample_sheet for organism etc. (sample-centric; same meta structure as runsheet)
        sample_sheet_rows = (is_tmt && params.sample_sheet)
            ? Channel.fromPath(params.sample_sheet).flatMap { samplesheetToList(it, "$projectDir/schema_sample_sheet.json") }
            : Channel.empty()

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

        if (params.uniprot_id && params.reference_proteome) {
            error "ERROR: Cannot specify both uniprot_id and reference_proteome."
        }

        proteome = Channel.empty()
        def staged_proteome = (ep == 'mzml') || params.uniprot_id || params.reference_proteome
        if (params.uniprot_id) {
            GET_PROTEOME(output_dir)
            proteome = GET_PROTEOME.out.proteome_fasta
            ch_published = ch_published.mix(pub(GET_PROTEOME.out.proteome_fasta, ch_root, 'Proteome'))
        } else if (params.reference_proteome) {
            CHECK_DECOYS_CONTAMS(output_dir, file(params.reference_proteome))
            proteome = CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked
            ch_published = ch_published.mix(pub(CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked, ch_root, 'Proteome'))
        } else if (ep == 'mzml' && params.reference_table) {
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
            ch_published = ch_published
                .mix(pub(FETCH_REFERENCE_PROTEOME.out.proteome_fasta, ch_root, 'Proteome'))
                .mix(pub(CHECK_DECOYS_CONTAMS.out.proteome_fasta_checked, ch_root, 'Proteome'))
        } else if (ep == 'mzml') {
            error "ERROR: Must specify uniprot_id, reference_proteome, or reference_table with a matching organism row."
        }

        // Generate manifest: TMT [data_sheet, sample_sheet]; LFQ [runsheet]
        ch_sample_sheet = (is_tmt && params.sample_sheet) ? Channel.fromPath(params.sample_sheet) : Channel.value(file("${projectDir}/bin/placeholder"))
        ch_sheets = sheet_fp.combine(ch_sample_sheet).map { s, ss -> [s, ss] }
        FRAGPIPE_METADATA_SETUP(output_dir, ch_sheets)
        manifest = FRAGPIPE_METADATA_SETUP.out.manifest
        ch_published = ch_published
            .mix(pub(FRAGPIPE_METADATA_SETUP.out.manifest, ch_root, 'Metadata'))
            .mix(pub(FRAGPIPE_METADATA_SETUP.out.experiment_annotation, ch_root, 'Metadata'))
            .mix(pub(FRAGPIPE_METADATA_SETUP.out.msstats_tmt_annotation, ch_root, 'Metadata'))
            .mix(pub(FRAGPIPE_METADATA_SETUP.out.sheets, ch_root, 'Metadata'))

        def fp_levels = params.fp_analyst_levels ?
            params.fp_analyst_levels.toString().split(',')*.trim().findAll { it } :
            (params.fragpipe_workflow?.startsWith('TMT') ? ['protein','gene','peptide','site'] : ['protein','peptide'])
        def tmt_quant = (params.fp_analyst_tmt_quant_type ?: 'abundance').toLowerCase().startsWith('ratio') ? 'ratio' : 'abundance'

        ch_msstats_csv = Channel.empty()
        ch_msstats_ptm = Channel.empty()
        ch_combined_protein = Channel.empty()
        ch_combined_peptide = Channel.empty()
        ch_combined_tables = Channel.empty()
        ch_tmt_report_tables = Channel.empty()
        ch_abundance_protein = Channel.empty()
        ch_abundance_peptide = Channel.empty()
        ch_abundance_gene = Channel.empty()
        ch_abundance_single_site = Channel.empty()
        ch_ratio_protein = Channel.empty()
        ch_ratio_peptide = Channel.empty()
        ch_ratio_gene = Channel.empty()
        ch_ratio_single_site = Channel.empty()
        ch_search_versions = Channel.empty()
        ch_pmultiqc_versions = Channel.empty()
        ch_vv = Channel.empty()

        if (ep == 'mzml') {
            // mzML only. MSCONVERT exists but is not wired — convert Thermo .raw upstream.
            samples_full.branch { meta ->
                remote: isRemoteSrc(meta)
                local: true
            }.set { ch_src }
            COPY_INPUT(ch_out_dir, ch_src.local.map { meta -> tuple(meta, file(meta.data_file.toString())) })
            FETCH_INPUT(ch_out_dir, ch_src.remote.map { meta -> tuple(meta, meta.data_file.toString()) })
            ch_mzml = COPY_INPUT.out.mzml_files.mix(FETCH_INPUT.out.mzml_files)
            RAWBEANS_QC_ALL(output_dir, ch_mzml.map { it[1] }.collect())
            ch_published = ch_published.mix(pub(RAWBEANS_QC_ALL.out.qc_report, ch_root, 'RawBeans'))

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

            FRAGPIPE_CONFIG_SETUP(output_dir, fragpipe_config, proteome)
            ch_published = ch_published.mix(pub(FRAGPIPE_CONFIG_SETUP.out.fragpipe_config, ch_root, 'Metadata'))
            // Join on id — do not combine(collect()) id lists: Groovy flattens List into tuple slots
            ch_fragpipe_mzml = keep_first_tr
                ? ch_mzml
                    .map { meta, mzml -> tuple(meta.id.toString(), mzml) }
                    .join(samples_fp.map { m -> tuple(m.id.toString(), 1) }, by: 0)
                    .map { id, mzml, fpMarker -> mzml }
                    .collect()
                : ch_mzml.map { it[1] }.collect()

            FRAGPIPE(output_dir, FRAGPIPE_CONFIG_SETUP.out.fragpipe_config, ch_fragpipe_tools, manifest, proteome, ch_fragpipe_mzml, FRAGPIPE_METADATA_SETUP.out.experiment_annotation)

            ///////////////////////////////////////////////////////////
            // END HEADLESS FRAGPIPE
            ///////////////////////////////////////////////////////////

            ch_msstats_csv = FRAGPIPE.out.msstats_csv
            ch_msstats_ptm = FRAGPIPE.out.msstats_ptm_csv
            ch_combined_protein = FRAGPIPE.out.combined_protein
            ch_combined_peptide = FRAGPIPE.out.combined_peptide
            ch_combined_tables = FRAGPIPE.out.combined_tables
            ch_tmt_report_tables = FRAGPIPE.out.tmt_report_tables
            ch_abundance_protein = FRAGPIPE.out.abundance_protein
            ch_abundance_peptide = FRAGPIPE.out.abundance_peptide
            ch_abundance_gene = FRAGPIPE.out.abundance_gene
            ch_abundance_single_site = FRAGPIPE.out.abundance_single_site
            ch_ratio_protein = FRAGPIPE.out.ratio_protein
            ch_ratio_peptide = FRAGPIPE.out.ratio_peptide
            ch_ratio_gene = FRAGPIPE.out.ratio_gene
            ch_ratio_single_site = FRAGPIPE.out.ratio_single_site
            ch_search_versions = COPY_INPUT.out.versions.mix(FETCH_INPUT.out.versions)
                .mix(RAWBEANS_QC_ALL.out.versions)
                .mix(FRAGPIPE.out.versions)

            ch_published = ch_published
                .mix(pub(FRAGPIPE.out.msstats_csv, ch_root, 'FragPipe'))
                .mix(pub(FRAGPIPE.out.msstats_ptm_csv, ch_root, 'FragPipe'))
            ch_vv = FRAGPIPE.out.msstats_csv.map { f -> tuple('fragpipe', f) }

            ch_fragpipe_output_dir = FRAGPIPE.out.fragpipe_manifest.map { fragpipe_manifest -> fragpipe_manifest.parent }
            PMULTIQC(output_dir.map { it + "/pmultiqc" }, ch_fragpipe_output_dir)
            ch_pmultiqc_versions = PMULTIQC.out.versions
            ch_published = ch_published
                .mix(pub(PMULTIQC.out.html, ch_root, 'pmultiqc'))
                .mix(pub(PMULTIQC.out.zipped_data, ch_root, 'pmultiqc'))

            ZIP_FRAGPIPE_OUTPUTS(ch_out_dir, ch_fragpipe_output_dir)
            ch_published = ch_published.mix(pub(ZIP_FRAGPIPE_OUTPUTS.out.zip, ch_root, 'FragPipe'))
        } else {
            if (!params.fragpipe_workflow) {
                error "entry_point=fragpipe_output requires --fragpipe_workflow (LFQ-MBR, TMT10, TMT16, TMT16-phospho)."
            }
            def fp_root = resolve_fp_root(params.fragpipe_output)
            def has_combined = existing_under(fp_root, ['combined_protein.tsv'])
            def has_abund = existing_under(fp_root, tmt_table_pats('protein', 'abundance'))
            if (is_tmt && !has_abund) {
                def extra = has_combined ? ' Folder looks LFQ (combined_protein.tsv).' : ''
                error "entry_point=fragpipe_output: no abundance_protein_*.tsv under ${fp_root}.${extra}"
            }
            if (!is_tmt && !has_combined) {
                def extra = has_abund ? ' Folder looks TMT (abundance_protein_*.tsv).' : ''
                error "entry_point=fragpipe_output: no combined_protein.tsv under ${fp_root}.${extra}"
            }
            if (is_tmt && tmt_quant == 'ratio') {
                require_under(fp_root, tmt_table_pats('protein', 'ratio'), 'ratio_protein_*.tsv')
            }
            if (fp_levels.contains('peptide')) {
                if (is_tmt) {
                    require_under(fp_root, tmt_table_pats('peptide', tmt_quant), "${tmt_quant}_peptide_*.tsv")
                } else {
                    require_under(fp_root, ['combined_peptide.tsv'], 'combined_peptide.tsv')
                }
            }
            if (is_tmt && fp_levels.contains('gene')) {
                require_under(fp_root, tmt_table_pats('gene', tmt_quant), "${tmt_quant}_gene_*.tsv")
            }
            if (is_tmt && fp_levels.contains('site')) {
                require_under(fp_root, tmt_table_pats('single-site', tmt_quant), "${tmt_quant}_single-site_*.tsv")
            }

            ch_msstats_csv = files_ch(existing_under(fp_root, ['msstats.csv']))
            ch_msstats_ptm = files_ch(existing_under(fp_root, ['msstats_ptm.csv']))
            ch_combined_protein = files_ch(existing_under(fp_root, ['combined_protein.tsv']))
            ch_combined_peptide = files_ch(existing_under(fp_root, ['combined_peptide.tsv']))
            ch_combined_tables = files_ch(existing_under(fp_root, ['combined*.tsv']))
            ch_tmt_report_tables = files_ch(existing_under(fp_root, ['tmt-report/*.tsv', 'abundance_*.tsv', 'ratio_*.tsv']))
            ch_abundance_protein = files_ch(existing_under(fp_root, tmt_table_pats('protein', 'abundance')))
            ch_abundance_peptide = files_ch(existing_under(fp_root, tmt_table_pats('peptide', 'abundance')))
            ch_abundance_gene = files_ch(existing_under(fp_root, tmt_table_pats('gene', 'abundance')))
            ch_abundance_single_site = files_ch(existing_under(fp_root, tmt_table_pats('single-site', 'abundance')))
            ch_ratio_protein = files_ch(existing_under(fp_root, tmt_table_pats('protein', 'ratio')))
            ch_ratio_peptide = files_ch(existing_under(fp_root, tmt_table_pats('peptide', 'ratio')))
            ch_ratio_gene = files_ch(existing_under(fp_root, tmt_table_pats('gene', 'ratio')))
            ch_ratio_single_site = files_ch(existing_under(fp_root, tmt_table_pats('single-site', 'ratio')))

            ch_vv = ch_msstats_csv.map { f -> tuple('fragpipe', f) }

            if (is_fp_workdir(fp_root)) {
                PMULTIQC(output_dir.map { it + "/pmultiqc" }, Channel.value(fp_root))
                ch_pmultiqc_versions = PMULTIQC.out.versions
                ch_published = ch_published
                    .mix(pub(PMULTIQC.out.html, ch_root, 'pmultiqc'))
                    .mix(pub(PMULTIQC.out.zipped_data, ch_root, 'pmultiqc'))
            }
        }

        // Pass in FragPipe tables to clean and publish; run downstream modules with original tables
        ch_fragpipe_tables_to_clean = (params.fragpipe_workflow?.startsWith('TMT') ?
            ch_tmt_report_tables :
            ch_combined_tables
        ).collect()
        CLEAN_FRAGPIPE_TABLES(
            ch_out_dir,
            ch_fragpipe_tables_to_clean,
            FRAGPIPE_METADATA_SETUP.out.experiment_annotation
        )
        ch_published = ch_published.mix(pub(CLEAN_FRAGPIPE_TABLES.out.tables, ch_root, 'FragPipe'))

        // Resolve gene_annotations_url for DE_results: direct file, or organism + reference table
        gene_annotations_url = Channel.value(null)
        if (params.gene_annotations_file) {
            gene_annotations_url = Channel.value(params.gene_annotations_file)
        } else if (params.reference_table) {
            gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
        }
        ch_lfq_versions = Channel.empty()
        ch_tmt_versions = Channel.empty()
        def do_drop = (params.drop_decoys_contams != false && params.drop_decoys_contams != 'false')
        def keep_fpar_contams = (params.fp_analyst_keep_contaminants == true || params.fp_analyst_keep_contaminants == 'true')

        // FRAGPIPEANALYSTR (FragPipeAnalystR): levels from params or workflow default. TMT uses tmt-report abundance/ratio; LFQ uses combined_*.
        ch_fp_analyst_inputs = Channel.empty()
        if (fp_levels.contains('protein')) {
            def ch_protein = params.fragpipe_workflow?.startsWith('TMT') ?
                (tmt_quant == 'ratio' ? ch_ratio_protein : ch_abundance_protein).map { takeFirst(it) } :
                ch_combined_protein
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_protein.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("protein", q, e) })
        }
        if (fp_levels.contains('peptide')) {
            def ch_peptide = params.fragpipe_workflow?.startsWith('TMT') ?
                (tmt_quant == 'ratio' ? ch_ratio_peptide : ch_abundance_peptide).map { takeFirst(it) } :
                ch_combined_peptide
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_peptide.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("peptide", q, e) })
        }
        if (params.fragpipe_workflow?.startsWith('TMT') && fp_levels.contains('gene')) {
            def ch_gene = (tmt_quant == 'ratio' ? ch_ratio_gene : ch_abundance_gene).map { takeFirst(it) }
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_gene.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("gene", q, e) })
        }
        if (params.fragpipe_workflow?.startsWith('TMT') && fp_levels.contains('site')) {
            def ch_site = (tmt_quant == 'ratio' ? ch_ratio_single_site : ch_abundance_single_site).map { takeFirst(it) }
            ch_fp_analyst_inputs = ch_fp_analyst_inputs.mix(
                ch_site.combine(FRAGPIPE_METADATA_SETUP.out.experiment_annotation).map { q, e -> tuple("site", q, e) })
        }

        ch_msstats_in = ch_msstats_csv
        if (do_drop) {
            ch_drop_in = ch_msstats_csv.map { f -> tuple('msstats', f) }
            if (!keep_fpar_contams) {
                ch_drop_in = ch_drop_in.mix(ch_fp_analyst_inputs.map { kind, q, e -> tuple(kind, q) })
            }
            DROP_DECOYS_CONTAMS(ch_drop_in)
            ch_msstats_in = DROP_DECOYS_CONTAMS.out.cleaned.filter { it[0] == 'msstats' }.map { it[1] }
            if (!keep_fpar_contams) {
                ch_fp_analyst_inputs = ch_fp_analyst_inputs
                    .map { kind, q, e -> tuple(kind, e) }
                    .join(DROP_DECOYS_CONTAMS.out.cleaned.filter { it[0] != 'msstats' })
                    .map { kind, e, q -> tuple(kind, q, e) }
            }
        }

        if (params.fragpipe_workflow == 'LFQ-MBR') {
            MSSTATS(output_dir, FRAGPIPE_METADATA_SETUP.out.experiment_annotation, ch_msstats_in)
            ch_lfq_versions = MSSTATS.out.versions
            ch_published = ch_published
                .mix(pub(MSSTATS.out.comparison, ch_root, 'MSstats'))
                .mix(pub(MSSTATS.out.contrasts, ch_root, 'MSstats'))
            ch_vv = ch_vv.mix(MSSTATS.out.comparison.map { f -> tuple('msstats', f) })
        }
        if (params.fragpipe_workflow?.startsWith('TMT')) {
            MSSTATSTMT(
                output_dir,
                FRAGPIPE_METADATA_SETUP.out.msstats_tmt_annotation,
                ch_msstats_in
            )
            ch_tmt_versions = MSSTATSTMT.out.versions
            ch_published = ch_published
                .mix(pub(MSSTATSTMT.out.comparison, ch_root, 'MSstatsTMT'))
                .mix(pub(MSSTATSTMT.out.contrasts, ch_root, 'MSstatsTMT'))
                .mix(pub(MSSTATSTMT.out.conditions_notice, ch_root, 'MSstatsTMT'))
                .mix(pub(MSSTATSTMT.out.runs_notice, ch_root, 'MSstatsTMT'))
            ch_vv = ch_vv.mix(MSSTATSTMT.out.comparison.map { f -> tuple('msstatstmt', f) })
        }

        ch_fp_with_annot = ch_fp_analyst_inputs.combine(gene_annotations_url)
        FRAGPIPEANALYSTR(ch_out_dir, ch_fp_with_annot)
        ch_published = ch_published.mix(pub_fpar(FRAGPIPEANALYSTR.out.published, ch_root))
        ch_vv = ch_vv.mix(
            FRAGPIPEANALYSTR.out.published.map { kind, files -> tuple("fpar_${kind}", files_only(as_list(files))) }
        )
        if (!params.skip_vv) {
            VV_STEP(ch_out_dir, ch_vv)
            VV_CONCAT_FILTER(ch_out_dir, VV_STEP.out.log.collect())
            ch_published = ch_published
                .mix(pub(VV_STEP.out.log, ch_root, 'VV_Logs'))
                .mix(pub(VV_CONCAT_FILTER.out.logs, ch_root, 'VV_Logs'))
        }

        ch_lfq_versions = ch_lfq_versions.mix(ch_tmt_versions).mix(FRAGPIPEANALYSTR.out.versions)

        // Software version capturing
        ch_software_versions = Channel.empty()
            | mix(FRAGPIPE_METADATA_SETUP.out.versions)
            | mix(ch_search_versions)
            | mix(ch_pmultiqc_versions)
            | mix(ch_lfq_versions)
            | mix(ch_tmt_versions)
        ch_software_versions
            | unique
            | collectFile(newLine: true)
            | set { ch_final_software_versions }
        SOFTWARE_VERSIONS(output_dir, ch_final_software_versions)

        if (params.uniprot_id) {
            ch_protocol_uniprot = Channel.value(params.uniprot_id.toString())
        } else if (staged_proteome && params.reference_table) {
            ch_protocol_uniprot = PARSE_ANNOTATIONS_TABLE.out.uniprot_id.map { uid -> uid?.toString() ?: '' }
        } else {
            ch_protocol_uniprot = Channel.value('')
        }
        ch_protocol_reference_table = Channel.value(staged_proteome ? (params.reference_table?.toString() ?: '') : '')
        ch_protocol_reference_proteome = Channel.value(staged_proteome ? (params.reference_proteome?.toString() ?: '') : '')
        ch_used_proteome = proteome.map { it.toString() }.ifEmpty('')
        GENERATE_PROCESSED_PROTOCOL(
            ch_out_dir,
            SOFTWARE_VERSIONS.out.software_versions,
            ch_used_proteome,
            ch_protocol_uniprot,
            ch_protocol_reference_table,
            ch_protocol_reference_proteome
        )
        ch_published = ch_published
            .mix(pub(SOFTWARE_VERSIONS.out.software_versions, ch_root, 'GeneLab'))
            .mix(pub(GENERATE_PROCESSED_PROTOCOL.out.processed_protocol, ch_root, 'GeneLab'))

    emit:
        published = ch_published
}
