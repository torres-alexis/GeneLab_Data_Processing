process GENERATE_METADATA_TABLE {
    // Generates a metadata table from an OSDR ISA archive and runsheet
    publishDir "${ ch_outdir }/Metadata",
        mode: params.publish_dir_mode,
        pattern: "*_metadata_table.txt"

    input:
        val(ch_outdir)
        val(glds_accession)
        path(isa_archive)
        path(runsheet)

    output:
        path("${ glds_accession }_metadata_table.txt"), emit: metadata_table
        path("versions.yaml"), emit: versions
    
    script:
        """
        create_table_v2.py --accession ${ glds_accession }  \
                        --isa-zip ${ isa_archive } \
                        --output-dir . \
                        --runsheet ${ runsheet }

        echo '"${task.process}":' > versions.yaml
        echo "    python: \$(python --version | sed 's/Python //')" >> versions.yaml
        """
}