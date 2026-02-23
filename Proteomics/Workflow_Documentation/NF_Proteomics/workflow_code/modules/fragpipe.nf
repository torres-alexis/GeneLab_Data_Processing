process FRAGPIPE {
    tag "${workflow_config.getName()}"
    containerOptions = "--cleanenv --bind \$PWD,\$HOME/.config,${projectDir}"
    
    publishDir "${output_dir}/FragPipe/",
        mode: params.publish_dir_mode,
        pattern: "output/**",
        saveAs: { filename -> filename.toString().replaceFirst(/^output\//, '') }
    publishDir "${output_dir}/Metadata/",
        mode: params.publish_dir_mode,
        pattern: "fragpipe-files.fp-manifest",
        saveAs: { filename -> "fragpipe-files.fp-manifest" }

    input:
    val(output_dir)
    path(workflow_config)
    path(fragpipe_tools)
    path(manifest)
    path(proteome)
    path(mzml_files)

    output:
    path("output/**"), emit: fragpipe_outputs
    path("fragpipe-files.fp-manifest"), emit: fragpipe_manifest
    path("output/msstats.csv"), emit: msstats_csv, optional: true
    path("output/experiment_annotation.tsv"), emit: experiment_annotation
    path("output/combined_protein.tsv"), emit: combined_protein_tsv
    path("output/combined_peptide.tsv"), emit: combined_peptide_tsv
    path("versions.yml"), emit: versions

    script:
    def workflow_config_basename = workflow_config.getName()
    def ram_gb = task.memory.toGiga().intValue()
    """
    # Export environment variables for FragPipe (as recommended in GitHub issue #755)
    export XDG_CONFIG_HOME=\${PWD}/fragpipe_home
    export JAVA_OPTS="-Djava.io.tmpdir=\${PWD}/fragpipe_temp"
    
    # Run FragPipe in headless mode
    /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe \\
        --headless \\
        --workflow ${workflow_config} \\
        --manifest ${manifest} \\
        --workdir . \\
        --ram ${ram_gb} \\
        --threads ${task.cpus} \\
        --config-tools-folder ${fragpipe_tools}
    #     --config-python /usr/bin/python3.11
    
    # After FragPipe runs, move everything from work folder into output folder
    # (except: mzML files, proteome, tools folder, manifest files, updated workflow config)
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
           [[ "\${item}" == fragpipe-files.fp-manifest ]] || \\
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
    FP_TOOLS=/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools
    TOOLS_DIR=${fragpipe_tools}
    LOG=\$(ls output/log_*.txt 2>/dev/null | head -1)
    echo '"'"${task.process}"'"': > versions.yml
    echo "    fragpipe: \$(/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe --help 2>&1 | grep -E '^FragPipe' | head -1 | sed 's|FragPipe v||')" >> versions.yml
    MSF_JAR=\$(find "\${TOOLS_DIR}" . -name 'MSFragger*.jar' 2>/dev/null | head -1)
    echo "    msfragger: \$(java -jar "\${MSF_JAR}" 2>&1 | grep 'MSFragger version' | cut -d' ' -f3 | sed 's|^MSFragger-||' || echo 'unknown')" >> versions.yml
    IQ_JAR=\$(find "\${TOOLS_DIR}" . -name 'IonQuant*.jar' 2>/dev/null | head -1)
    echo "    ionquant: \$([ -n "\${IQ_JAR}" ] && basename "\${IQ_JAR}" | sed 's|IonQuant-\\([0-9.]*\\)\\.jar|\\1|' || echo 'unknown')" >> versions.yml
    PHILO=\$(find "\${FP_TOOLS}" -path '*/Philosopher/philosopher-*' -type f 2>/dev/null | head -1)
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
