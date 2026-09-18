process ADD_DECOYS_CONTAMS {
    input:
    val(output_dir)
    path(proteome_fasta)

    output:
    path("output/*"), emit: proteome_fasta

    script:
    def decoy_tag = params.philosopher_decoy_prefix ?: 'rev_'
    def want_decoys = params.philosopher_decoys ? '--want-decoys' : ''
    def contam_prefix = (params.philosopher_contaminants && params.philosopher_contaminants_prefix) ? params.philosopher_contaminants_prefix.toString() : ''
    def decoy_flag = params.philosopher_decoys ? (decoy_tag != 'rev_' ? "--prefix ${decoy_tag}" : '') : '--nodecoys'
    def contam_flag = params.philosopher_contaminants ? (contam_prefix ? '--contam --contamprefix' : '--contam') : ''
    """
    set -euo pipefail
    input_fasta=${proteome_fasta}

    python3 \$(command -v add_decoys_contams.py) inspect "\$input_fasta" \\
        --decoy-prefix ${decoy_tag} \\
        --contam-prefix '${contam_prefix}' \\
        ${want_decoys} \\
        --emit-env inspect.env \\
        --emit-study-accs study_accs.txt
    . ./inspect.env

    mkdir -p output

    if [ "\$NEED_DECOY" != "true" ] && [ -z "${contam_flag}" ]; then
        echo "FASTA already has required decoys; copy through"
        cp "\$input_fasta" output/
        exit 0
    fi

    python3 \$(command -v add_decoys_contams.py) write-targets "\$input_fasta" \\
        --decoy-prefix ${decoy_tag} \\
        --output targets_only.fas

    FP_BASE=\$(ls -d /fragpipe_bin/fragpipe-*/fragpipe-*/ 2>/dev/null | head -1)
    [ -z "\$FP_BASE" ] && { echo "ERROR: FragPipe not found under /fragpipe_bin/fragpipe-*/" >&2; exit 1; }
    PHILO=\$(ls \${FP_BASE%/}/tools/Philosopher/philosopher-v* 2>/dev/null | head -1)
    [ -z "\$PHILO" ] && { echo "ERROR: Philosopher not found" >&2; exit 1; }

    "\${PHILO}" workspace --init
    "\${PHILO}" database --custom targets_only.fas ${decoy_flag} ${contam_flag}
    processed=\$(ls -1t -- *.fa* | grep -vxF -- targets_only.fas | head -1)
    [ -z "\$processed" ] && { echo "ERROR: Philosopher did not produce a FASTA output (*.fa*)" >&2; exit 1; }

    if [ -n "${contam_flag}" ] && [ -n "${contam_prefix}" ]; then
        python3 \$(command -v add_decoys_contams.py) restore-study "\$processed" \\
            --decoy-prefix ${decoy_tag} \\
            --contam-prefix '${contam_prefix}' \\
            --study-accs study_accs.txt \\
            --output output/"\$processed"
    else
        cp "\$processed" output/
    fi
    "\${PHILO}" workspace --clean
    """
}
