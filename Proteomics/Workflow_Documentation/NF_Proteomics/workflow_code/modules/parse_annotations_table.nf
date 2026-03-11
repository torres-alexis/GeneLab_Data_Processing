// Parse GeneLab Reference Annotations table to resolve organism -> gene_annotations_url.
// Table format: https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv
// Organism key: species column (index 1), e.g. "Homo sapiens". Input organism_sci: "homo_sapiens".

process PARSE_ANNOTATIONS_TABLE {
  tag "Organism: ${organism_sci}"

  input:
    val(annotations_csv_url_string)
    val(organism_sci)

  output:
    val(gene_annotations_url), emit: gene_annotations_url

  exec:
    def organisms = [:]
    if (annotations_csv_url_string.startsWith('http://') || annotations_csv_url_string.startsWith('https://')) {
      annotations_csv_url_string.toURL().splitEachLine(",") { fields ->
        organisms[fields[1]] = fields
      }
    } else {
      new File(annotations_csv_url_string).splitEachLine(",") { fields ->
        organisms[fields[1]] = fields
      }
    }

    def organism_key = organism_sci.capitalize().replace("_", " ")
    if (organisms.containsKey(organism_key)) {
      gene_annotations_url = organisms[organism_key][10]
      if (gene_annotations_url != null && gene_annotations_url.contains('figshare.com/ndownloader/files/')) {
        def file_id = (gene_annotations_url =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        gene_annotations_url = "https://api.figshare.com/v2/file/download/${file_id}"
      }
      println "Gene annotations URL for '${organism_key}': ${gene_annotations_url}"
    } else {
      println "WARNING: Organism '${organism_key}' not in annotations table. Skipping DE annotations."
      gene_annotations_url = null
    }
}
