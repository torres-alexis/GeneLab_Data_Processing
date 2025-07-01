{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/nextflow-io/nf-schema/refs/heads/master/examples/samplesheetToListMeta/pipeline/assets/schema_input.json",
  "title": "NF_Proteomics - params.runsheet schema",
  "description": "Schema for the file provided with params.runsheet",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "Sample Name": {
        "type": "string",
        "pattern": "^\\S+$",
        "errorMessage": "Sample name must be provided and cannot contain spaces",
        "meta": ["sample_id"]
      },
      "Raw Spectral Data File": {
        "type": "string",
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'",
        "meta": ["raw_data_file"]
      }
    },
    "required": ["Sample Name", "Raw Spectral Data File"]
  }
}