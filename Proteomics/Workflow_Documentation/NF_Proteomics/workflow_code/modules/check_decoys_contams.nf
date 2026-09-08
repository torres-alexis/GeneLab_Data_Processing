process CHECK_DECOYS_CONTAMS {
    publishDir path: { "${output_dir}/Proteome/" },
        mode: params.publish_dir_mode,
        pattern: "output/*",
        saveAs: { filename -> filename.replaceAll(/^output\//, '') }

    input:
    val(output_dir)
    path(proteome_fasta)

    output:
    path("output/*"), emit: proteome_fasta_checked

    script:
    def decoy_tag = params.philosopher_decoy_prefix ?: 'rev_'
    def needs_decoys = params.philosopher_decoys ? 'true' : 'false'
    def needs_contaminants = params.philosopher_contaminants ? 'true' : 'false'
    def contam_prefix = (params.philosopher_contaminants && params.philosopher_contaminants_prefix) ? 'true' : 'false'
    
    """
    input_fasta=${proteome_fasta}
    decoy_tag=${decoy_tag}
    
    # Check if decoys exist
    echo Checking decoys in \$input_fasta
    has_decoys=false
    if grep -q \$decoy_tag \$input_fasta; then
        has_decoys=true
        echo "Decoys found in the fasta file"
    fi
    
    echo Checking contaminants in \$input_fasta
    contaminant_ids="P02769|P00760|P00711|P13645|P04264|O43790|P00004|P00698|P01012|P02768|P99999|P32503"
    contaminant_pattern="\\|(\${contaminant_ids})\\|"
    has_contaminants=false
    if grep -qE "\$contaminant_pattern" \$input_fasta; then
        has_contaminants=true
        echo "Contaminants found in the fasta file"
    fi
    
    # Determine what needs to be added based on params
    needs_processing=false
    decoy_flag=""
    contam_flag=""
    
    if [ "${needs_decoys}" = "true" ] && [ "\$has_decoys" = "false" ]; then
        needs_processing=true
        decoy_flag="--prefix \$decoy_tag"
        echo "Decoys need to be added"
    fi
    
    if [ "${needs_contaminants}" = "true" ] && [ "\$has_contaminants" = "false" ]; then
        needs_processing=true
        if [ "${contam_prefix}" = "true" ]; then
            contam_flag="--contam --contamprefix"
        else
            contam_flag="--contam"
        fi
        echo "Contaminants need to be added"
    fi
    
    mkdir -p output
    
    if [ "\$needs_processing" = "true" ]; then
        echo "Processing file to add missing components"
        FP_BASE=\$(ls -d /fragpipe_bin/fragpipe-*/fragpipe-*/ 2>/dev/null | head -1)
        [ -z "\$FP_BASE" ] && { echo "ERROR: FragPipe not found under /fragpipe_bin/fragpipe-*/" >&2; exit 1; }
        PHILO=\$(ls \${FP_BASE%/}/tools/Philosopher/philosopher-v* 2>/dev/null | head -1)
        [ -z "\$PHILO" ] && { echo "ERROR: Philosopher not found" >&2; exit 1; }
        "\${PHILO}" workspace --init
        "\${PHILO}" database --custom \$input_fasta \$decoy_flag \$contam_flag
        processed_fasta=\$(ls -t *.fa* 2>/dev/null | head -1)
        [ -z "\$processed_fasta" ] && { echo "ERROR: Philosopher did not produce a FASTA output (*.fa*)" >&2; exit 1; }
        cp "\$processed_fasta" output/
        "\${PHILO}" workspace --clean
    else
        cp \$input_fasta output/
    fi
    """
}
