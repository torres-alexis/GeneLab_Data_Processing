// Parse GL-DPPD-7110-A reference annotations CSV for organism -> gene_annotations_url + proteome FASTA.
// Table: GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv
// Organism key: species column (index 1), e.g. "Mus musculus". Input organism_sci: "mus_musculus".

process PARSE_ANNOTATIONS_TABLE {
  tag "Organism: ${organism_sci}"

  input:
    val(annotations_csv_url_string)
    val(organism_sci)

  output:
    val(gene_annotations_url), emit: gene_annotations_url
    val(proteome_source), emit: proteome_source
    val(uniprot_id), emit: uniprot_id

  exec:
    def organisms = [:]
    if (annotations_csv_url_string.startsWith('http://') || annotations_csv_url_string.startsWith('https://')) {
      annotations_csv_url_string.toURL().splitEachLine(",") { fields ->
        if (fields.size() > 1) {
          organisms[fields[1]] = fields
        }
      }
    } else {
      new File(annotations_csv_url_string).splitEachLine(",") { fields ->
        if (fields.size() > 1) {
          organisms[fields[1]] = fields
        }
      }
    }

    def organism_key = organism_sci.capitalize().replace("_", " ")
    if (organisms.containsKey(organism_key)) {
      def row = organisms[organism_key]
      gene_annotations_url = row.size() > 10 ? row[10] : null
      uniprot_id = row.size() > 12 ? row[12] : null
      proteome_source = row.size() > 13 ? row[13] : null
      if (gene_annotations_url == '') gene_annotations_url = null
      if (uniprot_id == '') uniprot_id = null
      if (proteome_source == '') proteome_source = null
      println "Reference table match for '${organism_key}':"
      println "  gene_annotations_url: ${gene_annotations_url ?: '(empty)'}"
      println "  uniprot_id: ${uniprot_id ?: '(empty)'}"
      println "  proteome: ${proteome_source ?: '(empty)'}"
    } else {
      println "WARNING: Organism '${organism_key}' not in reference annotations table."
      gene_annotations_url = null
      uniprot_id = null
      proteome_source = null
    }
}
