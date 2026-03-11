process FRAGPIPE {
    tag "${workflow_config.getName()}"
    containerOptions = "--cleanenv --bind \$PWD,\$HOME/.config,${projectDir}"
    
    publishDir "${output_dir}/FragPipe/",
        mode: params.publish_dir_mode,
        pattern: "output/**",
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    // publishDir "${output_dir}/MSstats/",
    //     mode: params.publish_dir_mode,
    //     pattern: "output/msstats.csv",
    //     saveAs: { filename -> "msstats_fragpipe.csv" }
    // publishDir "${output_dir}/MSstats/",
    //     mode: params.publish_dir_mode,
    //     pattern: "output/msstats_ptm.csv",
    //     saveAs: { filename -> "msstats_ptm_fragpipe.csv" }

    input:
    val(output_dir)
    path(workflow_config)
    path(fragpipe_tools)
    path(manifest)
    path(proteome)
    path(mzml_files)
    path(experiment_annotation)

    output:
    path("output/**"), emit: fragpipe_outputs
    path("output/fragpipe-files.fp-manifest"), emit: fragpipe_manifest
    path("versions.yml"), emit: versions
    // LFQ / shared
    path("output/msstats.csv"), emit: msstats_csv, optional: true
    path("output/msstats_ptm.csv"), emit: msstats_ptm_csv, optional: true
    path("output/experiment_annotation.tsv"), emit: experiment_annotation, optional: true
    path("output/combined_protein.tsv"), emit: combined_protein, optional: true
    path("output/combined_peptide.tsv"), emit: combined_peptide, optional: true
    // TMT tmt-report/ (suffix from prot_norm: 0=None, 1=MD, 2=GN, -1=All; glob matches whatever TMT-Integrator produced)
    path("output/tmt-report/abundance_protein_*.tsv"), emit: abundance_protein, optional: true
    path("output/tmt-report/abundance_peptide_*.tsv"), emit: abundance_peptide, optional: true
    path("output/tmt-report/abundance_gene_*.tsv"), emit: abundance_gene, optional: true
    path("output/tmt-report/abundance_single-site_*.tsv"), emit: abundance_single_site, optional: true
    path("output/tmt-report/abundance_multi-site_*.tsv"), emit: abundance_multi_site, optional: true
    path("output/tmt-report/ratio_protein_*.tsv"), emit: ratio_protein, optional: true
    path("output/tmt-report/ratio_peptide_*.tsv"), emit: ratio_peptide, optional: true
    path("output/tmt-report/ratio_gene_*.tsv"), emit: ratio_gene, optional: true
    path("output/tmt-report/ratio_single-site_*.tsv"), emit: ratio_single_site, optional: true
    path("output/tmt-report/ratio_multi-site_*.tsv"), emit: ratio_multi_site, optional: true

    script:
    def workflow_config_basename = workflow_config.getName()
    def ram_gb = task.memory.toGiga().intValue()
    def tmt_label = params.fragpipe_workflow == 'TMT10' ? 'TMT10' :
        params.fragpipe_workflow == 'TMT16' || params.fragpipe_workflow == 'TMT16-phospho' ? 'TMT16' :
        params.fragpipe_workflow?.startsWith('TMT') ? 'TMT10' : 'TMT10'
    """
    # Export environment variables for FragPipe (as recommended in GitHub issue #755)
    export XDG_CONFIG_HOME=\${PWD}/fragpipe_home
    export JAVA_OPTS="-Djava.io.tmpdir=\${PWD}/fragpipe_temp"
    mkdir -p fragpipe_temp
    
    # TMT: reorganize mzML into plex-specific (manifest.tsv Experiment_Bioreplicate) folders + annotation.txt before FragPipe
    if [[ "${params.fragpipe_workflow}" == TMT* ]] && [[ -s ${experiment_annotation} ]]; then
        bash ${projectDir}/bin/tmt_stage_by_plex.sh ${manifest} ${experiment_annotation} ${tmt_label}
    fi
    
    FP_BASE=\$(ls -d /fragpipe_bin/fragpipe-*/fragpipe-*/ 2>/dev/null | head -1)
    [ -z "\$FP_BASE" ] && { echo "ERROR: FragPipe not found under /fragpipe_bin/fragpipe-*/" >&2; exit 1; }
    
    # Run FragPipe in headless mode
    "\${FP_BASE}bin/fragpipe" \\
        --headless \\
        --workflow ${workflow_config} \\
        --manifest ${manifest} \\
        --workdir . \\
        --ram ${ram_gb} \\
        --threads ${task.cpus} \\
        --config-tools-folder ${fragpipe_tools}
    
    # After FragPipe runs, move everything from work folder into output folder
    # (except: mzML files, proteome, tools folder, input manifest, updated workflow config)
    # Move entire folders/directories into output/, preserving structure
    mkdir output
    for item in *; do
        # Skip if it's one of the excluded files/folders
        if [[ "\${item}" == "output" ]] || \\
           [[ "\${item}" == "versions.yml" ]] || \\
           [[ "\${item}" == *.mzML ]] || \\
           [[ "\${item}" == *.fas ]] || \\
           [[ "\${item}" == "tools" ]] || \\
           [[ "\${item}" == manifest*.tsv ]] || \\
           [[ "\${item}" == "fragpipe_home" ]] || \\
           [[ "\${item}" == "fragpipe_temp" ]] || \\
           [[ "\${item}" == "${workflow_config_basename}" ]]; then
            continue
        fi
        # Move everything else (folders and files) to output
        if [[ -e "\${item}" ]]; then
            mv "\${item}" output/
        fi
    done
    
    # Export version info (FragPipe + bundled subtools). Use | as sed delimiter to avoid Groovy parsing //.
    FP_TOOLS="\${FP_BASE}tools"
    TOOLS_DIR=${fragpipe_tools}
    LOG=\$(ls output/log_*.txt 2>/dev/null | head -1)
    echo '"'"${task.process}"'"': > versions.yml
    echo "    fragpipe: \$("\${FP_BASE}bin/fragpipe" --help 2>&1 | grep -E '^FragPipe' | head -1 | sed 's|FragPipe v||')" >> versions.yml
    MSF_JAR=\$(find "\${TOOLS_DIR}" . -name 'MSFragger*.jar' 2>/dev/null | head -1)
    echo "    msfragger: \$(java -jar "\${MSF_JAR}" 2>&1 | grep 'MSFragger version' | cut -d' ' -f3 | sed 's|^MSFragger-||' || echo 'unknown')" >> versions.yml
    IQ_JAR=\$(find "\${TOOLS_DIR}" . -name 'IonQuant*.jar' 2>/dev/null | head -1)
    echo "    ionquant: \$([ -n "\${IQ_JAR}" ] && basename "\${IQ_JAR}" | sed 's|IonQuant-\\([0-9.]*\\)\\.jar|\\1|' || echo 'unknown')" >> versions.yml
    PHILO=\$(ls \${FP_TOOLS}/Philosopher/philosopher-v* 2>/dev/null | head -1)
    echo "    philosopher: \$([ -n "\${PHILO}" ] && "\${PHILO}" version 2>&1 | grep -oE 'version=v?[0-9.]+' | sed 's|version=v\\?||' || echo 'unknown')" >> versions.yml
    PERC=\$(find "\${FP_TOOLS}" -path '*/percolator_*/linux/percolator' -type f 2>/dev/null | head -1)
    echo "    percolator: \$([ -n "\${PERC}" ] && "\${PERC}" --help 2>&1 | grep -oE 'Percolator version [0-9.]+' | sed 's|Percolator version ||' || echo 'unknown')" >> versions.yml
    BATMASS_JAR=\$(find "\${FP_TOOLS}" -name 'batmass-io-*.jar' 2>/dev/null | head -1)
    echo "    batmass: \$([ -n "\${BATMASS_JAR}" ] && basename "\${BATMASS_JAR}" | sed 's|batmass-io-\\([0-9.]*\\)\\.jar|\\1|' || echo 'unknown')" >> versions.yml
    if [ -n "\${LOG}" ]; then
      echo "    msbooster: \$(grep 'MSBooster version' "\${LOG}" 2>/dev/null | head -1 | sed 's|MSBooster version ||' | tr -d '\\r' || echo 'unknown')" >> versions.yml
      echo "    diann: \$(grep 'DIA-NN version' "\${LOG}" 2>/dev/null | head -1 | sed 's|DIA-NN version ||' | tr -d '\\r' || echo 'unknown')" >> versions.yml
    fi
    """
}
