process STAGE_INPUT {
    tag "${meta.id}"
    
    publishDir "${output_dir}/00-RawData/${meta.id}",
        mode: params.publish_dir_mode,
        pattern: "*.mzML"

    input:
    val(output_dir)
    tuple val(meta), val(file_url)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml_files
    path("versions.yml"), emit: versions

    script:
    """
    # Stage file from various sources with retry logic
    echo "Staging ${meta.id} from: ${file_url}"
    
    # Detect and download/copy based on source type
    if [[ "${file_url}" =~ ^s3:// ]]; then
        echo "Detected S3 source"
        aws s3 cp "${file_url}" ./raw_file --retry-mode adaptive
    elif [[ "${file_url}" =~ ^https?:// ]]; then
        echo "Detected URL source"
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 -O raw_file "${file_url}"
    elif [[ -f "${file_url}" ]]; then
        echo "Detected local file"
        cp -P "${file_url}" ./raw_file
    else
        echo "ERROR: Unknown file source: ${file_url}"
        exit 1
    fi
    
    # Check if file was successfully staged
    if [[ ! -f raw_file ]]; then
        echo "ERROR: Failed to stage file from ${file_url}"
        exit 1
    fi
    
    # Detect format and handle accordingly
    file_type=\$(file raw_file)
    echo "File type detected: \$file_type"
    
    if [[ \$file_type == *"Zip archive"* ]]; then
        echo "Extracting ZIP archive..."
        unzip -q raw_file
        # Find mzML files in extracted content
        find . -name "*.mzML" -exec mv {} ${meta.id}${params.assay_suffix}.mzML \\;
    elif [[ \$file_type == *"gzip compressed"* ]]; then
        echo "Extracting gzipped file..."
        gunzip -c raw_file > ${meta.id}${params.assay_suffix}.mzML
    elif [[ "${file_url}" == *.mzML ]]; then
        echo "Already mzML format"
        mv raw_file ${meta.id}${params.assay_suffix}.mzML
    else
        echo "WARNING: Unknown format, assuming mzML"
        mv raw_file ${meta.id}${params.assay_suffix}.mzML
    fi
    
    # Verify final mzML file exists
    if [[ ! -f ${meta.id}${params.assay_suffix}.mzML ]]; then
        echo "ERROR: Failed to create standardized mzML file"
        exit 1
    fi
    
    echo "Successfully staged and standardized: ${meta.id}${params.assay_suffix}.mzML"
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    wget: \$(wget --version | head -n1 | cut -d' ' -f3)" >> versions.yml
    echo "    file: \$(file --version | head -n1)" >> versions.yml
    """
}
