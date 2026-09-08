process DROP_DECOYS_CONTAMS {
    tag "${kind}:${table.baseName}"

    input:
        tuple val(kind), path(table)

    output:
        tuple val(kind), path("cleaned/${table.name}"), emit: cleaned

    script:
    def decoy_prefix = (params.philosopher_decoy_prefix != null && params.philosopher_decoy_prefix != '') ? params.philosopher_decoy_prefix : 'rev_'
    """
    mkdir -p cleaned
    decoy_contam.py \\
        --input ${table} \\
        --output cleaned/${table.name} \\
        --decoy-prefix ${decoy_prefix}
    """
}
