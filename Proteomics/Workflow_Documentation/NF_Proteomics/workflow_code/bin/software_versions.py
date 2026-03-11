#!/usr/bin/env python
"""Consolidate process-emitted software versions into a markdown table and YAML.

Input: YAML file with format:
  "TaskName":
      "SoftwareName": version

Output: software_versions*.md (markdown table) and software_versions*.yaml
"""
from __future__ import annotations
import re
from pathlib import Path

import yaml
from packaging import version


CONFIG = {
    "proteomics": [
        ["NF_Proteomics", "https://github.com/nasa/GeneLab_Data_Processing/tree/master/Proteomics"],
        ["dp_tools", "https://github.com/J-81/dp_tools"],
        ["RawBeans", "https://bitbucket.org/incpm/prot-qc/src/master/protqc/"],
        ["Philosopher", "https://github.com/Nesvilab/philosopher/releases/latest"],
        ["FragPipe", "https://fragpipe.nesvilab.org/"],
        ["MultiQC", "https://multiqc.info/"],
        ["pmultiqc", "https://github.com/bigbio/pmultiqc"],
        ["R", "https://www.r-project.org/"],
        ["MSstats", "https://github.com/Vitek-Lab/MSstats"],
        ["FragPipeAnalystR", "https://github.com/Nesvilab/FragPipeAnalystR"],
    ]
}

# Skip these when processing (infra, not assay-specific)
SKIP_SOFTWARE = {"file", "wget", "python", "nextflow"}

# FragPipe-internal (bundled): exclude from software table
SKIP_FRAGPIPE_INTERNAL = {"batmass", "msfragger", "msbooster", "diann", "percolator", "ionquant"}


class NumericAsStringSafeLoader(yaml.SafeLoader):
    """YAML loader that does NOT coerce ints/floats; numeric-looking scalars stay as strings."""


for ch, resolvers in list(NumericAsStringSafeLoader.yaml_implicit_resolvers.items()):
    NumericAsStringSafeLoader.yaml_implicit_resolvers[ch] = [
        (tag, regexp)
        for tag, regexp in resolvers
        if tag not in ("tag:yaml.org,2002:int", "tag:yaml.org,2002:float")
    ]


def _represent_str_always_quoted(dumper, data: str):
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


class QuotedStringDumper(yaml.SafeDumper):
    """YAML dumper that always quotes strings to avoid downstream numeric coercion."""


QuotedStringDumper.add_representer(str, _represent_str_always_quoted)


def normalize_name(name: str, known_names: list) -> str:
    """Match software name against known names, ignoring case and special chars."""
    name_clean = re.sub(r"[^a-zA-Z0-9]", "", name.lower())
    for known, _ in known_names:
        if re.sub(r"[^a-zA-Z0-9]", "", known.lower()) == name_clean:
            return known
    print(f"Warning: Unknown software detected: {name}")
    return name


def compare_versions(v1, v2):
    """Safely compare version strings."""
    try:
        return version.parse(str(v1)) > version.parse(str(v2))
    except version.InvalidVersion:
        return str(v1) > str(v2)


def prefer_more_precise_representation(new_version: str, existing_version: str) -> bool:
    """When versions are equal numerically, prefer the one with more precision."""
    try:
        if version.parse(new_version) == version.parse(existing_version):
            new_digits = len(re.sub(r"\D", "", new_version))
            existing_digits = len(re.sub(r"\D", "", existing_version))
            if new_digits != existing_digits:
                return new_digits > existing_digits
            return len(new_version) > len(existing_version)
    except version.InvalidVersion:
        if new_version == existing_version:
            return False
        return len(new_version) > len(existing_version)
    return False


def main(
    versions_path: Path,
    output_path: Path,
    assay: str = "proteomics",
    workflow: str = None,
    workflow_version: str = None,
):
    software_urls = {name: url for name, url in CONFIG[assay]}
    known_names = CONFIG[assay]
    processed_versions = {}

    if workflow and workflow_version:
        processed_versions[workflow] = workflow_version

    with versions_path.open() as f:
        data = yaml.load(f, Loader=NumericAsStringSafeLoader)
        for task_info in data.values():
            if isinstance(task_info, dict):
                for software, ver in task_info.items():
                    if str(software).lower() in SKIP_SOFTWARE:
                        continue
                    normalized_name = normalize_name(software, known_names)
                    ver_str = str(ver).strip("'\"")
                    # RawBeans: container tag 1.6.4 ships protqc 1.6.3; report 1.6.4
                    if normalized_name == "RawBeans" and ver_str == "1.6.3":
                        ver_str = "1.6.4"
                    if normalized_name not in processed_versions:
                        processed_versions[normalized_name] = ver_str
                    else:
                        existing = processed_versions[normalized_name]
                        if compare_versions(ver_str, existing) or prefer_more_precise_representation(
                            ver_str, existing
                        ):
                            processed_versions[normalized_name] = ver_str

    if not processed_versions:
        print("No software versions found to process")
        return

    config_order = [name for name, _ in CONFIG[assay]]
    known_software = [x for x in config_order if x in processed_versions]
    unknown_software = sorted(
        [
            x
            for x in processed_versions
            if x not in config_order
            and re.sub(r"[^a-zA-Z0-9]", "", x.lower()) not in SKIP_FRAGPIPE_INTERNAL
        ]
    )
    ordered_programs = known_software + unknown_software

    # Build markdown table (no tabulate dependency)
    lines = ["| Program | Version | Relevant Links |", "| --- | --- | --- |"]
    for program in ordered_programs:
        ver = processed_versions[program]
        url = software_urls.get(program, "")
        lines.append(f"| {program} | {ver} | {url} |")
    output_path.write_text("\n".join(lines) + "\n")

    versions_dict = {p: processed_versions[p] for p in ordered_programs}
    yaml_output = output_path.with_suffix(".yaml")
    with yaml_output.open("w") as f:
        yaml.dump(versions_dict, f, sort_keys=False, default_flow_style=False, Dumper=QuotedStringDumper)
    print(f"Wrote {output_path} and {yaml_output}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Consolidate software versions into markdown table")
    parser.add_argument("input", type=Path, help="Path to combined versions YAML file")
    parser.add_argument("output", type=Path, help="Output markdown file path")
    parser.add_argument("--assay", choices=["proteomics"], default="proteomics")
    parser.add_argument("--workflow", type=str, help="Workflow name (e.g. NF_Proteomics)")
    parser.add_argument("--workflow_version", type=str, help="Workflow version")
    args = parser.parse_args()
    main(args.input, args.output, args.assay, args.workflow, args.workflow_version)
