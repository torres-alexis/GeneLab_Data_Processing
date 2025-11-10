#!/usr/bin/env python3
"""
Generate MSFragger config file from template and overrides.
"""
import sys
import re
from pathlib import Path

# Default MSFragger parameters (all set to None - will use template defaults)
DEFAULT_PARAMS = {
    'precursor_mass_lower': None,
    'precursor_mass_upper': None,
    'precursor_mass_units': None,
    'precursor_true_tolerance': None,
    'precursor_true_units': None,
    'precursor_mass_mode': None,
    'fragment_mass_tolerance': None,
    'fragment_mass_units': None,
    'calibrate_mass': None,
    'use_all_mods_in_first_search': None,
    'decoy_prefix': None,
    'deisotope': None,
    'deneutralloss': None,
    'isotope_error': None,
    'mass_offsets': None,
    'mass_offsets_detailed': None,
    'use_detailed_offsets': None,
    'remove_precursor_peak': None,
    'remove_precursor_range': None,
    'intensity_transform': None,
    'activation_types': None,
    'analyzer_types': None,
    'group_variable': None,
    'require_precursor': None,
    'reuse_dia_fragment_peaks': None,
    'write_calibrated_mzml': None,
    'write_uncalibrated_mzml': None,
    'write_mzbin_all': None,
    'mass_diff_to_variable_mod': None,
    'localize_delta_mass': None,
    'delta_mass_exclude_ranges': None,
    'fragment_ion_series': None,
    'ion_series_definitions': None,
    'labile_search_mode': None,
    'restrict_deltamass_to': None,
    'diagnostic_intensity_filter': None,
    'Y_type_masses': None,
    'diagnostic_fragments': None,
    'remainder_fragment_masses': None,
    'search_enzyme_name_1': None,
    'search_enzyme_cut_1': None,
    'search_enzyme_nocut_1': None,
    'search_enzyme_sense_1': None,
    'allowed_missed_cleavage_1': None,
    'search_enzyme_name_2': None,
    'search_enzyme_cut_2': None,
    'search_enzyme_nocut_2': None,
    'search_enzyme_sense_2': None,
    'allowed_missed_cleavage_2': None,
    'num_enzyme_termini': None,
    'clip_nTerm_M': None,
    'variable_mod_01': None,
    'variable_mod_02': None,
    'variable_mod_03': None,
    'variable_mod_04': None,
    'variable_mod_05': None,
    'variable_mod_06': None,
    'variable_mod_07': None,
    'variable_mod_08': None,
    'variable_mod_09': None,
    'variable_mod_10': None,
    'variable_mod_11': None,
    'variable_mod_12': None,
    'variable_mod_13': None,
    'variable_mod_14': None,
    'variable_mod_15': None,
    'variable_mod_16': None,
    'allow_multiple_variable_mods_on_residue': None,
    'max_variable_mods_per_peptide': None,
    'max_variable_mods_combinations': None,
    'output_format': None,
    'output_report_topN': None,
    'output_max_expect': None,
    'report_alternative_proteins': None,
    'precursor_charge': None,
    'override_charge': None,
    'digest_min_length': None,
    'digest_max_length': None,
    'digest_mass_range': None,
    'max_fragment_charge': None,
    'track_zero_topN': None,
    'zero_bin_accept_expect': None,
    'zero_bin_mult_expect': None,
    'check_spectral_files': None,
    'minimum_peaks': None,
    'use_topN_peaks': None,
    'min_fragments_modelling': None,
    'min_matched_fragments': None,
    'min_sequence_matches': None,
    'minimum_ratio': None,
    'clear_mz_range': None,
    'add_Cterm_peptide': None,
    'add_Nterm_peptide': None,
    'add_Cterm_protein': None,
    'add_Nterm_protein': None,
    'add_G_glycine': None,
    'add_A_alanine': None,
    'add_S_serine': None,
    'add_P_proline': None,
    'add_V_valine': None,
    'add_T_threonine': None,
    'add_C_cysteine': None,
    'add_L_leucine': None,
    'add_I_isoleucine': None,
    'add_N_asparagine': None,
    'add_D_aspartic_acid': None,
    'add_Q_glutamine': None,
    'add_K_lysine': None,
    'add_E_glutamic_acid': None,
    'add_M_methionine': None,
    'add_H_histidine': None,
    'add_F_phenylalanine': None,
    'add_R_arginine': None,
    'add_Y_tyrosine': None,
    'add_W_tryptophan': None,
    'add_B_user_amino_acid': None,
    'add_J_user_amino_acid': None,
    'add_O_user_amino_acid': None,
    'add_U_user_amino_acid': None,
    'add_X_user_amino_acid': None,
    'add_Z_user_amino_acid': None,
}

def parse_params_file(filepath):
    """Parse MSFragger params file into dict, preserving structure."""
    config = {}
    lines = []
    
    with open(filepath, 'r') as f:
        for line in f:
            lines.append(line.rstrip())
            # Skip comments and empty lines for parsing
            if not line.strip() or line.strip().startswith('#'):
                continue
            
            # Parse key = value format
            match = re.match(r'^([^=#]+?)\s*=\s*(.*?)(?:\s*#.*)?$', line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                config[key] = value
    
    return config, lines


def write_config(config, template_lines, output_file):
    """Write config file, preserving template structure and applying values."""
    output_lines = []
    
    for line in template_lines:
        # If it's a config line, replace the value
        if not line.strip() or line.strip().startswith('#'):
            output_lines.append(line)
        else:
            match = re.match(r'^([^=#]+?)\s*=\s*(.*?)(?:\s*(#.*))?$', line)
            if match:
                key = match.group(1).strip()
                comment = match.group(3) if match.group(3) else ''
                
                if key in config:
                    # Replace value, preserve comment
                    if comment:
                        output_lines.append(f"{key} = {config[key]}\t\t\t{comment}")
                    else:
                        output_lines.append(f"{key} = {config[key]}")
                else:
                    # Keep original line
                    output_lines.append(line)
            else:
                output_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines) + '\n')

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_msfragger_config.py <template_file> <output_file> [--key value ...]")
        sys.exit(1)
    
    template_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Start with default params (all None)
    overrides = DEFAULT_PARAMS.copy()
    
    # Parse command line arguments as --key value pairs
    i = 3
    while i < len(sys.argv):
        if sys.argv[i].startswith('--'):
            key = sys.argv[i][2:]  # Remove '--'
            if i + 1 < len(sys.argv):
                value = sys.argv[i + 1]
                if value and value != 'null':
                    overrides[key] = value
                i += 2
            else:
                i += 1
        else:
            i += 1
    
    # Parse template
    config, template_lines = parse_params_file(template_file)
    
    # Apply overrides (only non-None values override template)
    for key, value in overrides.items():
        if value is not None and value != 'null':
            config[key] = str(value)
    
    # Write output
    write_config(config, template_lines, output_file)

if __name__ == '__main__':
    main()

