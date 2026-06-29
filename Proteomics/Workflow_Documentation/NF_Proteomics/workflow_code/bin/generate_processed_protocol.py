#!/usr/bin/env python
"""Generate NF_Proteomics processed-data protocol text."""

import argparse
import json
import os
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="Generate NF_Proteomics processed-data protocol.")
    parser.add_argument("--from_protocol", default=None, help="Regenerate an existing protocol using its sidecar config JSON")
    parser.add_argument("--outdir", default=None, help="Output directory for protocol file")
    parser.add_argument("--software_table", default=None, help="Path to software versions markdown table")
    parser.add_argument("--assay_suffix", default="", help="e.g. _GLProteomics")
    parser.add_argument("--workflow_version", default="unknown", help="NF_Proteomics workflow version")
    parser.add_argument("--fragpipe_workflow", default="", help="FragPipe workflow preset")
    parser.add_argument("--fp_analyst_levels", default="", help="FragPipeAnalystR levels")
    parser.add_argument("--fp_analyst_zip", default="", help="Whether FragPipeAnalystR plot zips are enabled")
    parser.add_argument("--fp_analyst_lfq_type", default="", help="LFQ quantification type passed to FragPipeAnalystR")
    parser.add_argument("--fp_analyst_tmt_quant_type", default="", help="TMT quantification type passed to FragPipeAnalystR")
    parser.add_argument("--tmt_extraction_tool", default="", help="TMT reporter-ion extraction tool (Philosopher or IonQuant)")
    parser.add_argument("--normalization_method", default="", help="FragPipeAnalystR normalization method (stored in embedded config)")
    parser.add_argument("--imputation_type", default="", help="FragPipeAnalystR imputation method (stored in embedded config)")
    parser.add_argument("--de_alpha", default="", help="FragPipeAnalystR DE adjusted p-value threshold")
    parser.add_argument("--de_lfc", default="", help="FragPipeAnalystR DE minimum |log2 fold change|")
    parser.add_argument("--de_fdr", default="", help="FragPipeAnalystR multiple-testing correction method")
    parser.add_argument("--uniprot_id", default="", help="UniProt proteome ID, if used")
    parser.add_argument("--reference_proteome", default="", help="User-supplied reference proteome path param, if used")
    parser.add_argument("--reference_table", default="", help="GL-DPPD-7110-A reference table path/URL, if used for proteome")
    parser.add_argument("--used_proteome", default="", help="Proteome FASTA staged by the pipeline (decoys/contams applied)")
    parser.add_argument("--output", default=None, help="Protocol output filename")
    return parser.parse_args()


def parse_software_versions(md_file):
    versions = {}
    with open(md_file, "r") as handle:
        for line in handle:
            if not line.startswith("|") or line.startswith("| ---"):
                continue
            parts = [part.strip() for part in line.strip().strip("|").split("|")]
            if len(parts) >= 2 and parts[0] not in ("Program", ""):
                versions[parts[0]] = parts[1]
    return versions


def version(software_versions, program):
    return software_versions.get(program, "unknown")


def protocol_path_value(value):
    """Use basename for local paths; keep URLs unchanged for reproducible protocol metadata."""
    s = "" if value is None else str(value).strip()
    if not s or s.startswith(("http://", "https://")):
        return s
    return os.path.basename(s)


def protocol_document(fragpipe_workflow):
    docs = {
        "LFQ-MBR": "GL-DPPD-[LFQ-MBR]",
        "TMT10": "GL-DPPD-[TMT10]",
        "TMT16": "GL-DPPD-[TMT16]",
        "TMT16-phospho": "GL-DPPD-[TMT16-phospho]",
    }
    return docs.get(fragpipe_workflow, "GL-DPPD-[workflow]")


def preset_runs_msbooster(fragpipe_workflow):
    """GeneLab-checked FragPipe presets: msbooster.run-msbooster in conf/workflows/*.workflow."""
    if fragpipe_workflow in ("LFQ-MBR", "TMT10"):
        return True
    if fragpipe_workflow in ("TMT16", "TMT16-phospho"):
        return False
    return None


def default_fp_analyst_levels(fragpipe_workflow):
    if fragpipe_workflow == "LFQ-MBR":
        return "protein and peptide"
    if fragpipe_workflow == "TMT16-phospho":
        return "protein, gene, peptide, and site"
    if fragpipe_workflow in ("TMT10", "TMT16"):
        return "protein, gene, and peptide"
    return "configured"


def _str_or_default(val, default):
    s = "" if val is None else str(val).strip()
    return s if s else default


def tmt_extraction_tool_label(args):
    """GeneLab TMT preset uses Philosopher; matches params.tmt_extraction_tool / fragpipe_config_setup."""
    raw = _str_or_default(getattr(args, "tmt_extraction_tool", None), "Philosopher")
    return "IonQuant" if raw.strip().lower() == "ionquant" else "Philosopher"


def generate_protocol_content(args, software_versions):
    current_date = datetime.now().strftime("%Y-%m-%d")
    fragpipe_workflow = args.fragpipe_workflow or "configured workflow"
    fp_levels = args.fp_analyst_levels or default_fp_analyst_levels(fragpipe_workflow)
    dppd = protocol_document(fragpipe_workflow)
    fp_ver = version(software_versions, "FragPipe")
    de_alpha = _str_or_default(getattr(args, "de_alpha", None), "0.05")
    de_lfc = _str_or_default(getattr(args, "de_lfc", None), "1.0")
    de_fdr = _str_or_default(getattr(args, "de_fdr", None), "Benjamini Hochberg")

    header = f"# GeneLab Proteomics Pipeline Protocol — {dppd}{args.assay_suffix}\n"
    header += f"# Date: {current_date}\n\n"

    used_proteome = protocol_path_value(getattr(args, "used_proteome", None))
    used_proteome_note = f" ({used_proteome})" if used_proteome else ""

    database_sentence = ""
    if args.reference_table:
        uniprot_note = f" (UniProt proteome {args.uniprot_id})" if args.uniprot_id else ""
        database_sentence = (
            f"A pinned reference proteome from the GeneLab GL-DPPD-7110-A annotations table "
            f"({protocol_path_value(args.reference_table)}){uniprot_note}{used_proteome_note} was prepared "
            f"with decoys and contaminants using Philosopher (version {version(software_versions, 'Philosopher')}). "
        )
    elif args.uniprot_id:
        database_sentence = (
            f"A protein sequence database was generated from UniProt proteome {args.uniprot_id}"
            f"{used_proteome_note} with decoys and contaminants using Philosopher "
            f"(version {version(software_versions, 'Philosopher')}). "
        )
    elif args.reference_proteome:
        proteome_label = used_proteome or protocol_path_value(args.reference_proteome)
        database_sentence = (
            f"A user-supplied protein sequence database ({proteome_label}) was prepared "
            f"with decoys and contaminants using Philosopher (version {version(software_versions, 'Philosopher')}). "
        )
    elif used_proteome:
        database_sentence = (
            f"A protein sequence database ({used_proteome}) was prepared "
            f"with decoys and contaminants using Philosopher (version {version(software_versions, 'Philosopher')}). "
        )

    ptm_sentence = ""
    if fragpipe_workflow == "TMT16-phospho":
        ptm_sentence = "PTM site localization was performed with PTMProphet as part of the FragPipe workflow. "

    if fragpipe_workflow == "LFQ-MBR":
        quant_sentence = (
            f"Label-free MS1 quantification using {args.fp_analyst_lfq_type or 'Intensity'} values, "
            f"match-between-runs, and MaxLFQ summarization were performed with IonQuant. "
        )
        msstats_sentence = (
            f"Differential abundance analysis was performed with MSstats "
            f"(version {version(software_versions, 'MSstats')}). "
        )
    else:
        quant_type = args.fp_analyst_tmt_quant_type or "abundance"
        if fragpipe_workflow == "TMT16-phospho":
            level_list = "protein, peptide, gene, and modification-site"
        else:
            level_list = "protein, peptide, and gene"

        quant_sentence = (
            f"TMT reporter-ion intensities were extracted from tandem MS/MS spectra using {tmt_extraction_tool_label(args)}, "
            f"PSM-level quantification outputs were further processed with TMT-Integrator to generate summary "
            f"quantitative reports ({quant_type} quantification mode) at the {level_list} levels; integrated channel "
            f"abundances were log2 transformed and median centered. "
        )
        msstats_sentence = ""

    msb = preset_runs_msbooster(fragpipe_workflow)
    if msb is True:
        identification_sentence = (
            "Peptide-spectrum matching, deep learning-based spectral feature refinement, peptide-level rescoring, "
            "and protein inference were performed with MSFragger, MSBooster, Percolator, and Philosopher through "
            "FragPipe. "
        )
    else:
        identification_sentence = (
            "Peptide-spectrum matching, peptide-level rescoring, and protein inference were performed with "
            "MSFragger, Percolator, and Philosopher through FragPipe. "
        )

    description = (
        f"Data were processed as described in {dppd}, using NF_Proteomics version "
        f"{args.workflow_version}. In short, raw mass spectrometry files were staged as mzML files, "
        f"and raw data QC reports were generated with RawBeans (version {version(software_versions, 'RawBeans')}). "
        f"{database_sentence}"
        f"FragPipe (version {fp_ver}) was executed in headless (command-line) mode "
        f"with the \"{fragpipe_workflow}\" workflow preset. "
        f"{identification_sentence}"
        f"{ptm_sentence}"
        f"{quant_sentence}"
        f"QC metrics produced by FragPipe were summarized first using pmultiqc "
        f"(version {version(software_versions, 'pmultiqc')}); those summaries were then aggregated with MultiQC "
        f"(version {version(software_versions, 'MultiQC')}). "
        f"{msstats_sentence}"
        f"Downstream statistical analysis and visualizations were performed "
        f"at the {fp_levels} levels with FragPipeAnalystR (version {version(software_versions, 'FragPipeAnalystR')}), "
        f"including data quality control, limma-based differential expression analysis using adjusted p-value "
        f"threshold {de_alpha}, absolute log2 fold change greater than {de_lfc}, and {de_fdr} multiple-testing correction, "
        f"together with feature, volcano, and heatmap visualization, and gene ontology and pathway enrichment analysis."
    )

    return header + description + "\n"


def config_path_for_protocol(protocol_path):
    return f"{protocol_path}.config.json"


def embedded_config_from_protocol(protocol_path):
    start = "# NF_Proteomics protocol generation config:"
    lines = []
    in_config = False
    with open(protocol_path, "r") as handle:
        for line in handle:
            if line.rstrip("\n") == start:
                in_config = True
                continue
            if in_config:
                if not line.startswith("# "):
                    break
                lines.append(line[2:])
    if not lines:
        raise FileNotFoundError(config_path_for_protocol(protocol_path))
    return json.loads("".join(lines))


def load_saved_config(protocol_path):
    config_path = config_path_for_protocol(protocol_path)
    if os.path.exists(config_path):
        with open(config_path, "r") as handle:
            return json.load(handle)
    return embedded_config_from_protocol(protocol_path)


def save_config(args, output_path, software_versions):
    path_keys = {"reference_table", "reference_proteome", "used_proteome"}
    config = {
        key: protocol_path_value(value) if key in path_keys else value
        for key, value in vars(args).items()
        if key not in {"from_protocol", "outdir", "output", "software_table", "software_versions"}
    }
    config["software_versions"] = software_versions
    with open(config_path_for_protocol(output_path), "w") as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
    return config


def embedded_config_block(config):
    lines = ["\n# NF_Proteomics protocol generation config:\n"]
    for line in json.dumps(config, indent=2, sort_keys=True).splitlines():
        lines.append(f"# {line}\n")
    return "".join(lines)


def apply_saved_config(args):
    if not args.from_protocol:
        return args

    protocol_path = os.path.abspath(args.from_protocol)
    saved = load_saved_config(protocol_path)
    for key, value in saved.items():
        if getattr(args, key, None) in (None, ""):
            setattr(args, key, value)

    args.outdir = os.path.dirname(protocol_path)
    args.output = os.path.basename(protocol_path)
    return args


def main():
    args = apply_saved_config(parse_args())
    if not args.outdir:
        raise SystemExit("--outdir is required unless --from_protocol is used")
    if args.software_table:
        software_versions = parse_software_versions(args.software_table)
    else:
        software_versions = getattr(args, "software_versions", None)
        if not software_versions:
            raise SystemExit("--software_table is required unless --from_protocol uses embedded software_versions")

    os.makedirs(args.outdir, exist_ok=True)
    output = args.output or f"processed_data_protocol{args.assay_suffix}.txt"
    output_path = os.path.join(args.outdir, output)
    protocol_content = generate_protocol_content(args, software_versions)

    config = save_config(args, output_path, software_versions)
    with open(output_path, "w") as handle:
        handle.write(protocol_content)
        handle.write(embedded_config_block(config))

    print(f"Protocol file generated successfully: {output_path}")


if __name__ == "__main__":
    main()
