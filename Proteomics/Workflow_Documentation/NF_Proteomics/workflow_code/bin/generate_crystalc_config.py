#!/usr/bin/env python3
"""
Generate Crystal-C config file from template and overrides.
"""
import sys
import re
from pathlib import Path

def parse_args():
    """Parse command line arguments."""
    if len(sys.argv) < 3:
        print("Usage: generate_crystalc_config.py <template_file> <output_file> [--key value ...]", file=sys.stderr)
        sys.exit(1)
    
    template_file = sys.argv[1]
    output_file = sys.argv[2]
    overrides = {}
    
    # Parse --key value pairs
    i = 3
    while i < len(sys.argv):
        if sys.argv[i].startswith('--'):
            key = sys.argv[i][2:]  # Remove '--'
            if i + 1 < len(sys.argv):
                value = sys.argv[i + 1]
                overrides[key] = value
                i += 2
            else:
                print(f"Error: Missing value for {sys.argv[i]}", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"Error: Unexpected argument {sys.argv[i]}", file=sys.stderr)
            sys.exit(1)
    
    return template_file, output_file, overrides

def update_config_line(line, overrides):
    """Update a config line if it matches an override key."""
    # Match lines like: key = value  # comment
    match = re.match(r'^(\s*)(\w+)\s*=\s*(.+?)(\s*#.*)?$', line)
    if match:
        indent, key, value, comment = match.groups()
        if key in overrides:
            # Preserve comment if present
            comment_str = comment if comment else ''
            return f"{indent}{key} = {overrides[key]}{comment_str}\n"
    return line

def generate_config(template_file, output_file, overrides):
    """Generate config file from template with overrides."""
    template_path = Path(template_file)
    if not template_path.exists():
        print(f"Error: Template file not found: {template_file}", file=sys.stderr)
        sys.exit(1)
    
    with open(template_path, 'r') as f:
        template_lines = f.readlines()
    
    output_lines = []
    for line in template_lines:
        updated_line = update_config_line(line, overrides)
        output_lines.append(updated_line)
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.writelines(output_lines)

def main():
    template_file, output_file, overrides = parse_args()
    generate_config(template_file, output_file, overrides)

if __name__ == '__main__':
    main()

