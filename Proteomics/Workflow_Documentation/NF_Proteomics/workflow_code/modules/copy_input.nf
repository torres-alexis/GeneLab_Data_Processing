process COPY_INPUT {
    tag "${meta.id}"

    input:
    val(output_dir)
    tuple val(meta), path(src_file)

    output:
    tuple val(meta), path("staged/${meta.id}.mzML"), emit: mzml_files
    path("versions.yml"), emit: versions

    script:
    """
    echo "Copying ${meta.id} from: ${src_file}"
    mkdir -p staged
    cp -L "${src_file}" ./raw_file

    if unzip -t raw_file >/dev/null 2>&1; then
        unzip -q raw_file
        find . -name "*.mzML" ! -path "./staged/*" -exec mv {} staged/${meta.id}.mzML \\;
    elif gzip -t raw_file >/dev/null 2>&1; then
        gunzip -c raw_file > staged/${meta.id}.mzML
    else
        mv raw_file staged/${meta.id}.mzML
    fi

    if [[ ! -f staged/${meta.id}.mzML ]]; then
        echo "ERROR: Failed to create staged/${meta.id}.mzML"
        exit 1
    fi

    echo '"${task.process}":' > versions.yml
    """
}
