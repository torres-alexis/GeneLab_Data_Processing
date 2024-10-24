// modules/help.nf

def showHelp(workflow) {
    log.info """
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                         NASA GeneLab RNA-Seq Pipeline
                                Version: ${workflow.manifest.version}
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    Usage Examples:
    ----------------

    **Example 1:** Processing GLDS Datasets Using Genome FASTA and GTF from Ensembl
    ----------------------------------------------------------------------------
    ```
    nextflow run ./main.nf --osd 194 --glds 194
    ```

    **Example 2:** Processing GLDS Datasets Using Local Genome FASTA and GTF
    ---------------------------------------------------------------------
    ```
    nextflow run ./main.nf --osd 194 --glds 194 \\
        --reference_version 96 --reference_source <reference_label> \\
        --reference_fasta </path/to/fasta> --reference_gtf </path/to/gtf>
    ```

    **Example 3:** Processing Other Datasets with a User-Created Runsheet
    -----------------------------------------------------------------
    ```
    nextflow run ./main.nf --runsheet_path </path/to/runsheet>
    ```

    Parameters:
    ------------

    --help
        Show this help message and exit.

    --mode
        Workflow type (default: 'default'). Set to 'microbes' for processing microbes using Bowtie2.

    --osd
        OSD accession number to process (e.g., '194').

    --glds
        GLDS accession number to process (e.g., '194').

    --runsheet_path
        Path to a local runsheet instead of one automatically generated from a GLDS ISA archive.

    --isa_archive_path
        Path to a local ISA archive instead of retrieving from OSDR.

    --technical_replicates
        Path to a 2-column CSV file grouping duplicate samples. Example row: SampleA1,SampleA

    --truncateTo
        Subsample raw reads files to the specified number of reads for each raw reads file.

    --force_single_end
        Force analysis to use single-end processing. For paired-end datasets, only R1 is used. Default: false.

    --use_dummy_gene_counts
        Use random gene counts during DESeq2 (for testing purposes). Default: false.

    --reference_table
        URL or path to the reference table CSV file.

    --reference_store_path
        Directory where fetched reference files are downloaded. Default: './References'.

    --derived_store_path
        Directory where derived reference files are saved. Default: './DerivedReferences'.

    --reference_source
        Source of reference files (e.g., 'ensembl', 'ensembl_bacteria', 'ensembl_plants', 'ncbi').

    --reference_version
        Reference version (e.g., '96').

    --reference_fasta
        Path to a local reference FASTA file.

    --reference_gtf
        Path to a local reference GTF file.

    --validate_params
        Enable parameter validation. Default: true.

    --skip_vv
        Skip automated V&V. Default: false.

    --max_flag_code
        Maximum flag code. Default: 80.

    --rseqc_sample_count
        Number of reads to sample for RSeQC. Default: 15,000,000.

    --outdir
        Directory to save output files.

    --publish_dir_mode
        Publish directory mode. Default: 'link'.

    --email
        Email address for notifications.

    --version
        Show pipeline version and exit.

    --stage_local
        Whether to download the raw reads to the local filesystem. Default: true.

    Notes:
    -------

    - When specifying OSD and GLDS accession numbers, you only need to provide the numerical part (e.g., '194').
    - The pipeline will automatically prepend 'OSD-' and 'GLDS-' to these numbers internally.
    - Ensure that required parameters are provided for your specific use case.
    - For further assistance, consult the pipeline documentation or contact the development team.
    """.stripIndent()
    exit 0
}
