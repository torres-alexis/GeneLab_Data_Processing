def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(fasta_url), emit: reference_fasta_url
    val(gtf_url), emit: reference_gtf_url
    val(annotations_db_url), emit: annotations_db_url
    val(refSource), emit: reference_source
    val(refVersion), emit: reference_version
  
  exec:
    def organisms = [:]
    println "${colorCodes.c_yellow}Fetching table from ${annotations_csv_url_string}${colorCodes.c_reset}"
    
    // download data to memory
    annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    fasta_url = organisms[organism_key][5]
    gtf_url = organisms[organism_key][6]
    annotations_db_url = organisms[organism_key][9]
    refVersion = organisms[organism_key][3]
    refSource = organisms[organism_key][4]

    println "${colorCodes.c_blue}Annotation table values parsed for '${organism_key}':"
    println "${colorCodes.c_bright_green}- fasta_url: ${fasta_url}${colorCodes.c_reset}"
    println "${colorCodes.c_bright_green}- gtf_url: ${gtf_url}${colorCodes.c_reset}"
    println "${colorCodes.c_bright_green}- annotations_db_url: ${annotations_db_url}${colorCodes.c_reset}"
    println "${colorCodes.c_bright_green}- refVersion: ${refVersion}${colorCodes.c_reset}"
    println "${colorCodes.c_bright_green}- refSource: ${refSource}${colorCodes.c_reset}"
}
