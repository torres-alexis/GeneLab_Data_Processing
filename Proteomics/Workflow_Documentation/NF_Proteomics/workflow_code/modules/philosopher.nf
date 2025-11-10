process PHILOSOPHER {
    publishDir "${output_dir}/03-Philosopher/",
        mode: params.publish_dir_mode,
        pattern: "*.*"

    input:
    val(output_dir)
    path(database_fasta)
    path(pepxml_files)

    output:
    path("*.pep.xml"), emit: pepxml_validated
    path("*.prot.xml"), emit: protxml
    path("*.tsv"), emit: reports
    path("psm.tsv"), emit: psm_tsv
    path("protein.tsv"), emit: protein_tsv

    script:
    def fasta_file = database_fasta.toString()
    def decoy_prefix = params.philosopher_decoy_prefix ?: 'rev_'
    
    // Build peptideprophet command based on search type and Crystal-C status
    def peptideprophet_cmd
    if (params.search_type == "open") {
        if (params.run_crystal_c) {
            // Open search with Crystal-C (uses *_c.pepXML files)
            peptideprophet_cmd = "peptideprophet --nonparam --expectscore --decoyprobs --masswidth 1000.0 --clevel -2 --decoy ${decoy_prefix} --combine --database ${fasta_file} ./*_c.pepXML"
        } else {
            // Open search without Crystal-C (uses *.pepXML files)
            peptideprophet_cmd = "peptideprophet --nonparam --expectscore --decoyprobs --masswidth 1000.0 --clevel -2 --decoy ${decoy_prefix} --combine --database ${fasta_file} ./*.pepXML"
        }
    } else if (params.search_type == "nonspecific") {
        // Non-specific closed search
        peptideprophet_cmd = "peptideprophet --nonparam --expectscore --decoyprobs --ppm --accmass --nontt --decoy ${decoy_prefix} --database ${fasta_file} ./*.pepXML"
    } else {
        // Closed search (default)
        peptideprophet_cmd = "peptideprophet --nonparam --expectscore --decoyprobs --ppm --accmass --decoy ${decoy_prefix} --database ${fasta_file} ./*.pepXML"
    }
    
    // Build filter command based on search type
    def filter_cmd
    if (params.search_type == "open") {
        // Open search uses interact.pep.xml
        filter_cmd = "filter --sequential --razor --mapmods --tag ${decoy_prefix} --pepxml ./interact.pep.xml --protxml ./combined.prot.xml"
    } else {
        // Closed or non-specific closed search
        filter_cmd = "filter --sequential --razor --mapmods --tag ${decoy_prefix} --pepxml ./ --protxml ./combined.prot.xml"
    }
    
    """
    # Run PeptideProphet, ProteinProphet, and FDR filtering with Philosopher
    # See: https://fragpipe.nesvilab.org/docs/tutorial_headless.html
    
    # Philosopher is at /fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2
    philosopher_path="/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2"
    
    # Clean and initialize workspace
    \${philosopher_path} workspace --clean
    \${philosopher_path} workspace --init
    
    # Annotate database (sets up workspace metadata)
    \${philosopher_path} database --annotate ${fasta_file} --prefix ${decoy_prefix}
    
    # Run PeptideProphet (command varies by search type and Crystal-C status)
    # PeptideProphet automatically handles multiple files correctly:
    # - With --combine: processes all files together (single model)
    # - Without --combine: processes each file separately (separate models)
    \${philosopher_path} ${peptideprophet_cmd}
    
    # Run ProteinProphet
    # For open searches with --combine, PeptideProphet produces interact.pep.xml (singular)
    # For closed/nonspecific searches, PeptideProphet produces *.pep.xml files
    if [ "${params.search_type}" == "open" ]; then
        \${philosopher_path} proteinprophet --maxppmdiff 2000000 --output combined ./interact.pep.xml
    else
        \${philosopher_path} proteinprophet --maxppmdiff 2000000 --output combined ./*.pep.xml
    fi
    
    # Run Filter (command varies by search type)
    \${philosopher_path} ${filter_cmd}
    
    # Generate reports
    \${philosopher_path} report
    
    # Clean workspace
    \${philosopher_path} workspace --clean
    """
}

