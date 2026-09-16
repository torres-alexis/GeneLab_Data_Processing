process GENERATE_MD5SUMS {
    tag { type }

    input:
        tuple val(ch_outdir), val(type), path(processing_info)

    output:
        path("*md5sum${params.assay_suffix}.tsv"), emit: md5sum

    script:
        def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
        """
        generate_md5sums.py --outdir ${ch_outdir} ${assay_suffix_arg} --${type}
        """
}
