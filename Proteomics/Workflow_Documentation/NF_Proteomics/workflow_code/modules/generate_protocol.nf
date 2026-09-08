process GENERATE_PROCESSED_PROTOCOL {
    input:
        path(ch_outdir)
        path(software_versions_md)
        path(proteome_fasta)
        val(uniprot_id)
        val(reference_table)
        val(reference_proteome)

    output:
        path("processed_data_protocol${params.assay_suffix}.txt"), emit: processed_protocol

    script:
        """
        generate_processed_protocol.py \\
            --outdir . \\
            --software_table ${software_versions_md} \\
            --assay_suffix "${params.assay_suffix}" \\
            --workflow_version "${workflow.manifest.version}" \\
            --fragpipe_workflow "${params.fragpipe_workflow ?: ''}" \\
            --fp_analyst_levels "${params.fp_analyst_levels ?: ''}" \\
            --fp_analyst_zip "${params.fp_analyst_zip ?: ''}" \\
            --fp_analyst_lfq_type "${params.fp_analyst_lfq_type ?: ''}" \\
            --fp_analyst_tmt_quant_type "${params.fp_analyst_tmt_quant_type ?: ''}" \\
            --tmt_extraction_tool "${params.tmt_extraction_tool ?: 'Philosopher'}" \\
            --normalization_method "${params.fp_analyst_normalization_method ?: ''}" \\
            --imputation_type "${params.fp_analyst_imputation_type ?: ''}" \\
            --de_alpha "${params.fp_analyst_de_alpha != null ? params.fp_analyst_de_alpha : ''}" \\
            --de_lfc "${params.fp_analyst_de_lfc != null ? params.fp_analyst_de_lfc : ''}" \\
            --de_fdr "${params.fp_analyst_de_fdr ?: ''}" \\
            --uniprot_id "${uniprot_id}" \\
            --reference_proteome "${reference_proteome}" \\
            --reference_table "${reference_table}" \\
            --used_proteome ${proteome_fasta} \\
            --min_appearance_one_condition "${params.fp_analyst_min_appearance_one_cond ?: 50}" \\
            --min_global_appearance "${params.fp_analyst_min_global_appearance ?: 0}" \\
            --drop_decoys_contams "${params.drop_decoys_contams != false && params.drop_decoys_contams != 'false' ? 'true' : 'false'}" \\
            --output "processed_data_protocol${params.assay_suffix}.txt"
        """
}
