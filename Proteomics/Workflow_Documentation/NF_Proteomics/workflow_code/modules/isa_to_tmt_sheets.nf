// STUB: ISA to TMT data_sheet + sample_sheet. To be implemented in dp_tools (dpt-isa-to-tmt-sheets).
// For now errors. Use --data_sheet and --sample_sheet for TMT.
process ISA_TO_TMT_SHEETS {
    tag "${osd_accession}"

    input:
    val(ch_outdir)
    val(osd_accession)
    val(glds_accession)
    path(isa_archive)
    path(dp_tools_plugin)

    output:
    path("data_sheet.csv"), emit: data_sheet
    path("sample_sheet.csv"), emit: sample_sheet

    script:
    """
    echo "ERROR: ISA to TMT sheets not yet implemented. Use --data_sheet and --sample_sheet." >&2
    echo "Future: dpt-isa-to-tmt-sheets --accession \${glds_accession} --isa-archive \${isa_archive} --plugin-dir \${dp_tools_plugin}" >&2
    exit 1
    """
}
