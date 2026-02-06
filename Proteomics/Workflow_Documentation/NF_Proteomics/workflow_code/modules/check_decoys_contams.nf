process CHECK_DECOYS_CONTAMS {
    tag "${proteome_fasta.baseName}"
    publishDir "${output_dir}/Proteome/",
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
    
    # Check if contaminants exist (look for common contaminant UniProt IDs from cRAP)
    # Selected IDs: Laboratory (P02769, P00760, P00711), Dust/Contact (P13645, P04264, O43790), 
    # MW markers (P00004, P00698, P01012), UPS (P02768, P99999), Viral (P32503)
    echo Checking contaminants in \$input_fasta
    contaminant_ids="P02769|P00760|P00711|P13645|P04264|O43790|P00004|P00698|P01012|P02768|P99999|P32503"
    contaminant_pattern="\\|(${contaminant_ids})\\|"
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
        contam_flag="--contam"
        echo "Contaminants need to be added"
    fi
    
    mkdir -p output
    
    if [ "\$needs_processing" = "true" ]; then
        echo "Processing file to add missing components"
        /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 workspace --init
        /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 database --custom \$input_fasta \$decoy_flag \$contam_flag
        /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 workspace --clean
    fi
    
    cp \$input_fasta output/
    """
}
