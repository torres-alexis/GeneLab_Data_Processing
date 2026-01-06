#!/usr/bin/env python3

"""
Setup FragPipe workflow config file (.workflow):
- Add/update database.db-path with proteome path
- Apply assay suffix to output filename
- Convert to JSON format
"""

import argparse
import json
import sys
import os
from pathlib import Path


def setup_workflow_config(input_file, proteome_path, assay_suffix, output_dir):
    """Read workflow config, add/update database.db-path, and save with assay suffix."""
    # Read the input workflow file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Get proteome basename (just the filename)
    proteome_basename = os.path.basename(proteome_path)
    
    # Prepare the database.db-path line
    db_path_line = f"database.db-path={proteome_basename}\n"
    
    # Find if database.db-path already exists
    db_path_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith('database.db-path='):
            db_path_index = i
            break
    
    # Add or update database.db-path
    if db_path_index is not None:
        # Update existing entry
        lines[db_path_index] = db_path_line
    else:
        # Add after the first comment line (usually line 0 or 1)
        # Find first non-comment, non-empty line
        insert_index = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('#'):
                insert_index = i + 1
            elif line.strip():
                # Found first non-comment line, insert before it
                insert_index = i
                break
        lines.insert(insert_index, db_path_line)
    
    # Get input basename and apply assay suffix
    input_basename = os.path.basename(input_file)
    input_name, input_ext = os.path.splitext(input_basename)
    
    if assay_suffix:
        output_basename = f"{input_name}{assay_suffix}{input_ext}"
    else:
        output_basename = input_basename
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write the modified workflow file
    output_path = os.path.join(output_dir, output_basename)
    with open(output_path, 'w') as f:
        f.writelines(lines)
    
    return output_path, lines


def workflow_to_json(lines):
    """Convert workflow config lines to JSON dict."""
    config = {}
    
    for line in lines:
        # Skip comments and empty lines
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        # Parse key=value
        if '=' in line:
            parts = line.split('=', 1)  # Split only on first =
            key = parts[0].strip()
            value_str = parts[1].strip() if len(parts) > 1 else ''
            
            if key:
                config[key] = value_str
    
    return config


def main():
    parser = argparse.ArgumentParser(
        description="Setup FragPipe workflow config: add proteome path and convert to JSON"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input workflow config file (.workflow)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for modified workflow file"
    )
    parser.add_argument(
        "--proteome",
        required=True,
        help="Path to proteome FASTA file (basename will be used for database.db-path)"
    )
    parser.add_argument(
        "--assay_suffix",
        default="",
        help="Suffix to append to output filename (e.g., _GLProteomics)"
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="Print JSON to stdout"
    )
    args = parser.parse_args()
    
    try:
        # Setup workflow config (add proteome path, apply assay suffix)
        output_path, modified_lines = setup_workflow_config(
            args.input,
            args.proteome,
            args.assay_suffix,
            args.output
        )
        
        # Convert to JSON
        config = workflow_to_json(modified_lines)
        
        # Output JSON to stdout if requested
        if args.stdout:
            print(json.dumps(config, indent=2))
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

