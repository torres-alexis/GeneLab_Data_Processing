// Stage pinned reference proteome FASTA from GL-DPPD-7110-A annotations table (proteome column).
// Figshare ndownloader URLs rewritten like RNAseq annotate_dge_table.R → api.figshare.com/v2/file/download/{id}
// Download pattern matches RNAseq download_references.nf (wget keeps remote filename).
process FETCH_REFERENCE_PROTEOME {

    input:
    val(output_dir)
    val(proteome_source)

    output:
    path("*.fa*"), emit: proteome_fasta

    script:
    def src = proteome_source.toString()
    def download_url = src
    def figshare = (src =~ /figshare\.com\/ndownloader\/files\/([0-9]+)/)
    if (figshare.find()) {
        download_url = "https://api.figshare.com/v2/file/download/${figshare.group(1)}"
    }
    """
    src='${src}'
    download_url='${download_url}'
    echo "Staging reference proteome from: \$src"
    if [[ "\$download_url" != "\$src" ]]; then
        echo "Figshare ndownloader URL rewritten to: \$download_url"
    fi

    mkdir -p temp_fasta

    if [[ "\$src" =~ ^s3:// ]]; then
        echo "Detected S3 source"
        aws s3 cp "\$src" "temp_fasta/\$(basename "\$src")" --retry-mode adaptive
    elif [[ "\$download_url" =~ ^https?:// ]]; then
        echo "Detected URL source"
        wget --content-disposition --trust-server-names \\
            --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \\
            --directory-prefix temp_fasta "\$download_url"
    elif [[ -f "\$src" ]]; then
        echo "Detected local file"
        cp -P "\$src" "temp_fasta/\$(basename "\$src")"
    else
        echo "ERROR: Unknown proteome source: \$src"
        exit 1
    fi

    if ! ls temp_fasta/* &> /dev/null; then
        echo "ERROR: Failed to stage proteome from \$src"
        exit 1
    fi

    if ls temp_fasta/*.gz &> /dev/null; then
        echo "Extracting gzipped FASTA..."
        gunzip temp_fasta/*.gz
    fi

    mv temp_fasta/* ./

    if ! ls *.fa* &> /dev/null; then
        echo "ERROR: Failed to stage proteome FASTA"
        exit 1
    fi

    rm -rf temp_fasta
    echo "Successfully staged: \$(ls *.fa* | head -1)"
    """
}
