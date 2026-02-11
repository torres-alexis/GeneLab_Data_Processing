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
        --config-tools-folder ${fragpipe_tools} \\
        --config-python /usr/bin/python3.11 
    #     --ram 0 \\
    #     --threads -1
    
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
    
    # Export version info
    echo '"${task.process}":' > versions.yml
    echo "    fragpipe: \$(/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe --help 2>&1 | grep -E '^FragPipe' | head -1 | sed 's/FragPipe v//')" >> versions.yml
    """
}
