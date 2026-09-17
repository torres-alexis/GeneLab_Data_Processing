process FETCH_INPUT {
    tag "${meta.id}"

    input:
    val(output_dir)
    tuple val(meta), val(src)

    output:
    tuple val(meta), path("staged/${meta.id}.mzML"), emit: mzml_files
    path("versions.yml"), emit: versions

    script:
    """
    echo "Fetching ${meta.id} from: ${src}"
    mkdir -p staged

    if [[ "${src}" =~ ^s3:// ]]; then
        aws s3 cp "${src}" ./raw_file --retry-mode adaptive
    elif [[ "${src}" =~ ^https?:// ]]; then
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 -O raw_file "${src}"
    else
        echo "ERROR: FETCH_INPUT requires s3:// or http(s):// (${src})"
        exit 1
    fi

    if [[ ! -f raw_file ]]; then
        echo "ERROR: Failed to fetch ${src}"
        exit 1
    fi

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
