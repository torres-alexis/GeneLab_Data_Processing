process STAGE_INPUT {
    tag "${meta.id}"

    input:
    val(output_dir)
    val(meta)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml_files
    path("versions.yml"), emit: versions

    script:
    """
    echo "Staging ${meta.id} from: ${meta.data_file}"

    if [[ "${meta.data_file}" =~ ^s3:// ]]; then
        echo "Detected S3 source"
        aws s3 cp "${meta.data_file}" ./raw_file --retry-mode adaptive
    elif [[ "${meta.data_file}" =~ ^https?:// ]]; then
        echo "Detected URL source"
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 -O raw_file "${meta.data_file}"
    elif [[ -f "${meta.data_file}" ]]; then
        echo "Detected local file"
        cp -P "${meta.data_file}" ./raw_file
    else
        echo "ERROR: Unknown file source: ${meta.data_file}"
        exit 1
    fi
    
    if [[ ! -f raw_file ]]; then
        echo "ERROR: Failed to stage file from ${meta.data_file}"
        exit 1
    fi
    
    if unzip -t raw_file >/dev/null 2>&1; then
        echo "Extracting ZIP archive..."
        unzip -q raw_file
        find . -name "*.mzML" -exec mv {} ${meta.id}.mzML \\;
    elif gzip -t raw_file >/dev/null 2>&1; then
        echo "Extracting gzipped file..."
        gunzip -c raw_file > ${meta.id}.mzML
    elif [[ "${meta.data_file}" == *.mzML ]]; then
        echo "Already mzML format"
        mv raw_file ${meta.id}.mzML
    else
        echo "WARNING: Unknown format, assuming mzML"
        mv raw_file ${meta.id}.mzML
    fi
    
    if [[ ! -f ${meta.id}.mzML ]]; then
        echo "ERROR: Failed to create standardized mzML file"
        exit 1
    fi
    
    echo "Successfully staged and standardized: ${meta.id}.mzML"
    
    echo '"${task.process}":' > versions.yml
    echo "    wget: \$(wget --version | head -n1 | cut -d' ' -f3)" >> versions.yml
    """
}
