process GENERATE_MD5SUMS {
    input:
        val(ch_outdir)
        path(processing_info)

    output:
        path("raw_md5sum${params.assay_suffix}.tsv"), emit: raw_md5sum, optional: true
        path("processed_md5sum${params.assay_suffix}.tsv"), emit: processed_md5sum

    script:
        """
        ${projectDir}/bin/generate_md5sums.py --outdir ${ch_outdir} --assay_suffix ${params.assay_suffix}
        """
}