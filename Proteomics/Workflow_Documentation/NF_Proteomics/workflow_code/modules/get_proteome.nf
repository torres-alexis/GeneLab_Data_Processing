process GET_PROTEOME {
    tag "${params.uniprot_id}"
    publishDir "${output_dir}/Proteome/",
        mode: params.publish_dir_mode,
        pattern: "*.fas"

    input:
    val(output_dir)

    output:
    path("*.fas"), emit: proteome_fasta

    script:
    // Build flags based on params (accessing directly from params since they're globals)
    def reviewed_flag = params.philosopher_reviewed ? '--reviewed' : ''
    def isoforms_flag = params.philosopher_isoforms ? '--isoform' : ''
    def contaminants_flag = params.philosopher_contaminants ? '--contam' : ''
    def contaminants_prefix_flag = (params.philosopher_contaminants && params.philosopher_contaminants_prefix) ? '--contamprefix' : ''
    def nodecoys_flag = params.philosopher_decoys ? '' : '--nodecoys'
    def decoy_prefix_flag = params.philosopher_decoy_prefix && params.philosopher_decoy_prefix != 'rev_' ? "--prefix ${params.philosopher_decoy_prefix}" : ''
    def enzyme_flag = params.philosopher_enzyme && params.philosopher_enzyme != 'trypsin' ? "--enzyme ${params.philosopher_enzyme}" : ''
    def spike_in_flag = params.philosopher_spike_in ? "--add ${params.philosopher_spike_in}" : ''
    
    """
    # Use Philosopher to download and prepare proteome database from UniProt
    # See: https://github.com/Nesvilab/philosopher/wiki/Database
    
    # Philosopher is at /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2
    # Clean and initialize workspace (required by Philosopher)
    /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 workspace --clean
    /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 workspace --init
    
    # Download proteome using philosopher database command with --id
    /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2 database \\
        --id ${params.uniprot_id} \\
        ${reviewed_flag} \\
        ${isoforms_flag} \\
        ${contaminants_flag} \\
        ${contaminants_prefix_flag} \\
        ${decoy_prefix_flag} \\
        ${nodecoys_flag} \\
        ${enzyme_flag} \\
        ${spike_in_flag}
    """
}

