process IONQUANT {
    publishDir "${output_dir}/04-IonQuant/",
        mode: params.publish_dir_mode,
        pattern: "*.*"

    input:
    val(output_dir)
    path(psm_tsv_file)
    path(protein_tsv_file)
    path(mzml_files)

    output:
    path("*.tsv"), emit: quantification
    path("*.csv"), emit: msstats
    path("versions.yml"), emit: versions

    script:
    def memory_gb = task.memory.toGiga()
    def multidir_arg = params.ionquant_multidir ? "--multidir ${params.ionquant_multidir}" : (params.ionquant_msstats == 1 ? "--multidir ." : "")
    """
    # Run IonQuant for protein/peptide quantification
    java -Xmx${memory_gb}G -jar /usr/local/share/ionquant-1.11.9-0/IonQuant.jar \
        --specdir . \
        --psm ${psm_tsv_file} \
        --perform-ms1quant ${params.ionquant_perform_ms1quant} \
        --perform-isoquant ${params.ionquant_perform_isoquant} \
        --isotol ${params.ionquant_isotol} \
        --isolevel ${params.ionquant_isolevel} \
        --isotype ${params.ionquant_isotype} \
        ${params.ionquant_annotation ? "--annotation ${params.ionquant_annotation}" : ""} \
        --site-reports ${params.ionquant_site_reports} \
        --msstats ${params.ionquant_msstats} \
        --threads ${task.cpus} \
        --mztol ${params.ionquant_mztol} \
        --imtol ${params.ionquant_imtol} \
        --rttol ${params.ionquant_rttol} \
        ${multidir_arg} \
        --normalization ${params.ionquant_normalization} \
        --minisotopes ${params.ionquant_minisotopes} \
        --minscans ${params.ionquant_minscans} \
        --minions ${params.ionquant_minions} \
        ${params.ionquant_excludemods ? "--excludemods ${params.ionquant_excludemods}" : ""} \
        --maxlfq ${params.ionquant_maxlfq} \
        --ibaq ${params.ionquant_ibaq} \
        --minexps ${params.ionquant_minexps} \
        --minfreq ${params.ionquant_minfreq} \
        --tp ${params.ionquant_tp} \
        --mbr ${params.ionquant_mbr} \
        --mbrrttol ${params.ionquant_mbrrttol} \
        --mbrimtol ${params.ionquant_mbrimtol} \
        --mbrtoprun ${params.ionquant_mbrtoprun} \
        --mbrmincorr ${params.ionquant_mbrmincorr} \
        --ionmobility ${params.ionquant_ionmobility} \
        --ionfdr ${params.ionquant_ionfdr} \
        --peptidefdr ${params.ionquant_peptidefdr} \
        --proteinfdr ${params.ionquant_proteinfdr} \
        ${params.ionquant_light ? "--light ${params.ionquant_light}" : ""} \
        ${params.ionquant_medium ? "--medium ${params.ionquant_medium}" : ""} \
        ${params.ionquant_heavy ? "--heavy ${params.ionquant_heavy}" : ""} \
        ${params.ionquant_formula ? "--formula ${params.ionquant_formula}" : ""} \
        --requantify ${params.ionquant_requantify} \
        --writeindex ${params.ionquant_writeindex} \
        --locprob ${params.ionquant_locprob} \
        ${params.ionquant_modlist ? "--modlist ${params.ionquant_modlist}" : ""} \
        --uniqueness ${params.ionquant_uniqueness} \
        --intensitymode ${params.ionquant_intensitymode}
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    ionquant: \$(java -jar /usr/local/share/ionquant-1.11.9-0/IonQuant.jar --version 2>&1 | grep 'IonQuant version' | sed 's/IonQuant version //' || echo 'unknown')" >> versions.yml
    echo "    java: \$(java -version 2>&1 | head -n1 | cut -d' ' -f3 | tr -d '\"')" >> versions.yml
    """
}

