process MSSTATSTMTPTM {
    input:
    val(output_dir)
    path(msstats_tmt_annotation)
    path(msstats_csv)
    path(fragpipe_workflow_config)
    path(msstats_protein_csv, stageAs: 'msstats_protein.csv')
    path(msstats_protein_annotation, stageAs: 'msstats_protein_annotation.csv')

    output:
    path("versions.yml"), emit: versions
    path("msstatstmtptm_comparison*.csv"), emit: comparison, optional: true
    path("msstatstmtptm_ptm_comparison*.csv"), emit: ptm_comparison, optional: true
    path("msstatstmtptm_protein_comparison*.csv"), emit: protein_comparison, optional: true
    path("msstatstmtptm_adjusted_comparison*.csv"), emit: adjusted_comparison, optional: true
    path("msstatstmtptm_contrasts*.csv"), emit: contrasts, optional: true
    path("dropped-conditions-msstatstmtptm*.txt"), emit: conditions_notice, optional: true
    path("dropped-runs-msstatstmtptm*.txt"), emit: runs_notice, optional: true

    script:
    """
    sed -i '1s/,Probability,/,PeptideProphetProbability,/' ${msstats_csv}

    prot_csv=""
    prot_annot=""
    if [ -f msstats_protein.csv ] && [ -s msstats_protein.csv ]; then
        sed -i '1s/,Probability,/,PeptideProphetProbability,/' msstats_protein.csv
        prot_csv=msstats_protein.csv
    fi
    if [ -f msstats_protein_annotation.csv ] && [ -s msstats_protein_annotation.csv ]; then
        prot_annot=msstats_protein_annotation.csv
    fi

    fragpipe_ptm_mod_id.py ${fragpipe_workflow_config} > mod_ids.txt
    while IFS= read -r mod_id; do
        [ -z "\${mod_id}" ] && continue
        if [ -n "\${prot_csv}" ]; then
            msstatstmtptm_analysis.R . ${msstats_tmt_annotation} ${msstats_csv} ${params.assay_suffix} "\${mod_id}" "\${prot_csv}" "\${prot_annot}"
        else
            msstatstmtptm_analysis.R . ${msstats_tmt_annotation} ${msstats_csv} ${params.assay_suffix} "\${mod_id}"
        fi
    done < mod_ids.txt

    echo '"${task.process}":' > versions.yml
    echo "    MSstatsPTM: \$(Rscript -e 'cat(as.character(packageVersion(\"MSstatsPTM\")))' 2>/dev/null || echo 'unknown')" >> versions.yml
    echo "    r: \$(R --version 2>&1 | head -n1 | sed 's/.*version \\([0-9.]*\\).*/\\1/' || echo 'unknown')" >> versions.yml
    """
}
