process MSFRAGGER {
    publishDir "${output_dir}/01-MSFragger/",
        mode: params.publish_dir_mode,
        pattern: "*.{pepXML,pin,pepindex}"

    input:
    val(output_dir)
    path(mzml_file)
    path(msfragger_config)
    path(proteome)

    output:
    path("*.pepXML"), emit: pepxml
    path("*.pin"), emit: pin
    path("*.pepindex"), emit: pepindex
    path("versions.yml"), emit: versions

    script:
    def memory_gb = task.memory.toGiga()
    """
    
    # Define search type: open, closed, or nonspecific
    search_type="${params.search_type ?: 'closed'}"

    # MSFragger run for batch of mzML files
    java -Xmx${memory_gb}G -jar /usr/local/share/msfragger-4.2-0/MSFragger.jar \\
        ${msfragger_config} \\
        ${mzml_file}
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    msfragger: \$(java -jar /usr/local/share/msfragger-4.2-0/MSFragger.jar 2>&1 | grep 'MSFragger version' | cut -d' ' -f3)" >> versions.yml
    echo "    java: \$(java -version 2>&1 | head -n1 | cut -d' ' -f3 | tr -d '\"')" >> versions.yml
    """
}

process MSFRAGGER_CONFIG {
    publishDir "${output_dir}/01-MSFragger/",
        mode: params.publish_dir_mode,
        pattern: "msfragger_config.params"

    input:
    val(output_dir)
    val(search_type)
    path(template_file)
    path(database_fasta)

    output:
    path("msfragger_config.params"), emit: msfragger_config

    script:
    def database_filename = database_fasta.getName()
    
    // Only override database_name and num_threads (database_name should be just filename, not path)
    def override_args = " --database_name '${database_filename}' --num_threads ${task.cpus}"
    if (params.msfragger_precursor_mass_units != null) override_args += " --precursor_mass_units '${params.msfragger_precursor_mass_units}'"
    if (params.msfragger_precursor_true_tolerance != null) override_args += " --precursor_true_tolerance '${params.msfragger_precursor_true_tolerance}'"
    if (params.msfragger_precursor_true_units != null) override_args += " --precursor_true_units '${params.msfragger_precursor_true_units}'"
    if (params.msfragger_fragment_mass_tolerance != null) override_args += " --fragment_mass_tolerance '${params.msfragger_fragment_mass_tolerance}'"
    if (params.msfragger_fragment_mass_units != null) override_args += " --fragment_mass_units '${params.msfragger_fragment_mass_units}'"
    if (params.msfragger_calibrate_mass != null) override_args += " --calibrate_mass '${params.msfragger_calibrate_mass}'"
    if (params.msfragger_use_all_mods_in_first_search != null) override_args += " --use_all_mods_in_first_search '${params.msfragger_use_all_mods_in_first_search}'"
    if (params.msfragger_decoy_prefix != null) override_args += " --decoy_prefix '${params.msfragger_decoy_prefix}'"
    if (params.msfragger_deisotope != null) override_args += " --deisotope '${params.msfragger_deisotope}'"
    if (params.msfragger_deneutralloss != null) override_args += " --deneutralloss '${params.msfragger_deneutralloss}'"
    if (params.msfragger_isotope_error != null) override_args += " --isotope_error '${params.msfragger_isotope_error}'"
    if (params.msfragger_mass_offsets != null) override_args += " --mass_offsets '${params.msfragger_mass_offsets}'"
    if (params.msfragger_mass_offsets_detailed != null) override_args += " --mass_offsets_detailed '${params.msfragger_mass_offsets_detailed}'"
    if (params.msfragger_use_detailed_offsets != null) override_args += " --use_detailed_offsets '${params.msfragger_use_detailed_offsets}'"
    if (params.msfragger_precursor_mass_mode != null) override_args += " --precursor_mass_mode '${params.msfragger_precursor_mass_mode}'"
    if (params.msfragger_remove_precursor_peak != null) override_args += " --remove_precursor_peak '${params.msfragger_remove_precursor_peak}'"
    if (params.msfragger_remove_precursor_range != null) override_args += " --remove_precursor_range '${params.msfragger_remove_precursor_range}'"
    if (params.msfragger_intensity_transform != null) override_args += " --intensity_transform '${params.msfragger_intensity_transform}'"
    if (params.msfragger_activation_types != null) override_args += " --activation_types '${params.msfragger_activation_types}'"
    if (params.msfragger_analyzer_types != null) override_args += " --analyzer_types '${params.msfragger_analyzer_types}'"
    if (params.msfragger_group_variable != null) override_args += " --group_variable '${params.msfragger_group_variable}'"
    if (params.msfragger_require_precursor != null) override_args += " --require_precursor '${params.msfragger_require_precursor}'"
    if (params.msfragger_reuse_dia_fragment_peaks != null) override_args += " --reuse_dia_fragment_peaks '${params.msfragger_reuse_dia_fragment_peaks}'"
    if (params.msfragger_write_calibrated_mzml != null) override_args += " --write_calibrated_mzml '${params.msfragger_write_calibrated_mzml}'"
    if (params.msfragger_write_uncalibrated_mzml != null) override_args += " --write_uncalibrated_mzml '${params.msfragger_write_uncalibrated_mzml}'"
    if (params.msfragger_write_mzbin_all != null) override_args += " --write_mzbin_all '${params.msfragger_write_mzbin_all}'"
    if (params.msfragger_mass_diff_to_variable_mod != null) override_args += " --mass_diff_to_variable_mod '${params.msfragger_mass_diff_to_variable_mod}'"
    if (params.msfragger_localize_delta_mass != null) override_args += " --localize_delta_mass '${params.msfragger_localize_delta_mass}'"
    if (params.msfragger_delta_mass_exclude_ranges != null) override_args += " --delta_mass_exclude_ranges '${params.msfragger_delta_mass_exclude_ranges}'"
    if (params.msfragger_fragment_ion_series != null) override_args += " --fragment_ion_series '${params.msfragger_fragment_ion_series}'"
    if (params.msfragger_ion_series_definitions != null) override_args += " --ion_series_definitions '${params.msfragger_ion_series_definitions}'"
    if (params.msfragger_labile_search_mode != null) override_args += " --labile_search_mode '${params.msfragger_labile_search_mode}'"
    if (params.msfragger_restrict_deltamass_to != null) override_args += " --restrict_deltamass_to '${params.msfragger_restrict_deltamass_to}'"
    if (params.msfragger_diagnostic_intensity_filter != null) override_args += " --diagnostic_intensity_filter '${params.msfragger_diagnostic_intensity_filter}'"
    if (params.msfragger_Y_type_masses != null) override_args += " --Y_type_masses '${params.msfragger_Y_type_masses}'"
    if (params.msfragger_diagnostic_fragments != null) override_args += " --diagnostic_fragments '${params.msfragger_diagnostic_fragments}'"
    if (params.msfragger_remainder_fragment_masses != null) override_args += " --remainder_fragment_masses '${params.msfragger_remainder_fragment_masses}'"
    if (params.msfragger_search_enzyme_name_1 != null) override_args += " --search_enzyme_name_1 '${params.msfragger_search_enzyme_name_1}'"
    if (params.msfragger_search_enzyme_cut_1 != null) override_args += " --search_enzyme_cut_1 '${params.msfragger_search_enzyme_cut_1}'"
    if (params.msfragger_search_enzyme_nocut_1 != null) override_args += " --search_enzyme_nocut_1 '${params.msfragger_search_enzyme_nocut_1}'"
    if (params.msfragger_search_enzyme_sense_1 != null) override_args += " --search_enzyme_sense_1 '${params.msfragger_search_enzyme_sense_1}'"
    if (params.msfragger_allowed_missed_cleavage_1 != null) override_args += " --allowed_missed_cleavage_1 '${params.msfragger_allowed_missed_cleavage_1}'"
    if (params.msfragger_search_enzyme_name_2 != null) override_args += " --search_enzyme_name_2 '${params.msfragger_search_enzyme_name_2}'"
    if (params.msfragger_search_enzyme_cut_2 != null) override_args += " --search_enzyme_cut_2 '${params.msfragger_search_enzyme_cut_2}'"
    if (params.msfragger_search_enzyme_nocut_2 != null) override_args += " --search_enzyme_nocut_2 '${params.msfragger_search_enzyme_nocut_2}'"
    if (params.msfragger_search_enzyme_sense_2 != null) override_args += " --search_enzyme_sense_2 '${params.msfragger_search_enzyme_sense_2}'"
    if (params.msfragger_allowed_missed_cleavage_2 != null) override_args += " --allowed_missed_cleavage_2 '${params.msfragger_allowed_missed_cleavage_2}'"
    if (params.msfragger_num_enzyme_termini != null) override_args += " --num_enzyme_termini '${params.msfragger_num_enzyme_termini}'"
    if (params.msfragger_clip_nTerm_M != null) override_args += " --clip_nTerm_M '${params.msfragger_clip_nTerm_M}'"
    if (params.msfragger_variable_mod_01 != null) override_args += " --variable_mod_01 '${params.msfragger_variable_mod_01}'"
    if (params.msfragger_variable_mod_02 != null) override_args += " --variable_mod_02 '${params.msfragger_variable_mod_02}'"
    if (params.msfragger_variable_mod_03 != null) override_args += " --variable_mod_03 '${params.msfragger_variable_mod_03}'"
    if (params.msfragger_variable_mod_04 != null) override_args += " --variable_mod_04 '${params.msfragger_variable_mod_04}'"
    if (params.msfragger_variable_mod_05 != null) override_args += " --variable_mod_05 '${params.msfragger_variable_mod_05}'"
    if (params.msfragger_variable_mod_06 != null) override_args += " --variable_mod_06 '${params.msfragger_variable_mod_06}'"
    if (params.msfragger_variable_mod_07 != null) override_args += " --variable_mod_07 '${params.msfragger_variable_mod_07}'"
    if (params.msfragger_variable_mod_08 != null) override_args += " --variable_mod_08 '${params.msfragger_variable_mod_08}'"
    if (params.msfragger_variable_mod_09 != null) override_args += " --variable_mod_09 '${params.msfragger_variable_mod_09}'"
    if (params.msfragger_variable_mod_10 != null) override_args += " --variable_mod_10 '${params.msfragger_variable_mod_10}'"
    if (params.msfragger_variable_mod_11 != null) override_args += " --variable_mod_11 '${params.msfragger_variable_mod_11}'"
    if (params.msfragger_variable_mod_12 != null) override_args += " --variable_mod_12 '${params.msfragger_variable_mod_12}'"
    if (params.msfragger_variable_mod_13 != null) override_args += " --variable_mod_13 '${params.msfragger_variable_mod_13}'"
    if (params.msfragger_variable_mod_14 != null) override_args += " --variable_mod_14 '${params.msfragger_variable_mod_14}'"
    if (params.msfragger_variable_mod_15 != null) override_args += " --variable_mod_15 '${params.msfragger_variable_mod_15}'"
    if (params.msfragger_variable_mod_16 != null) override_args += " --variable_mod_16 '${params.msfragger_variable_mod_16}'"
    if (params.msfragger_allow_multiple_variable_mods_on_residue != null) override_args += " --allow_multiple_variable_mods_on_residue '${params.msfragger_allow_multiple_variable_mods_on_residue}'"
    if (params.msfragger_max_variable_mods_per_peptide != null) override_args += " --max_variable_mods_per_peptide '${params.msfragger_max_variable_mods_per_peptide}'"
    if (params.msfragger_max_variable_mods_combinations != null) override_args += " --max_variable_mods_combinations '${params.msfragger_max_variable_mods_combinations}'"
    if (params.msfragger_output_format != null) override_args += " --output_format '${params.msfragger_output_format}'"
    if (params.msfragger_output_report_topN != null) override_args += " --output_report_topN '${params.msfragger_output_report_topN}'"
    if (params.msfragger_output_max_expect != null) override_args += " --output_max_expect '${params.msfragger_output_max_expect}'"
    if (params.msfragger_report_alternative_proteins != null) override_args += " --report_alternative_proteins '${params.msfragger_report_alternative_proteins}'"
    if (params.msfragger_precursor_charge != null) override_args += " --precursor_charge '${params.msfragger_precursor_charge}'"
    if (params.msfragger_override_charge != null) override_args += " --override_charge '${params.msfragger_override_charge}'"
    if (params.msfragger_digest_min_length != null) override_args += " --digest_min_length '${params.msfragger_digest_min_length}'"
    if (params.msfragger_digest_max_length != null) override_args += " --digest_max_length '${params.msfragger_digest_max_length}'"
    if (params.msfragger_digest_mass_range != null) override_args += " --digest_mass_range '${params.msfragger_digest_mass_range}'"
    if (params.msfragger_max_fragment_charge != null) override_args += " --max_fragment_charge '${params.msfragger_max_fragment_charge}'"
    if (params.msfragger_track_zero_topN != null) override_args += " --track_zero_topN '${params.msfragger_track_zero_topN}'"
    if (params.msfragger_zero_bin_accept_expect != null) override_args += " --zero_bin_accept_expect '${params.msfragger_zero_bin_accept_expect}'"
    if (params.msfragger_zero_bin_mult_expect != null) override_args += " --zero_bin_mult_expect '${params.msfragger_zero_bin_mult_expect}'"
    if (params.msfragger_check_spectral_files != null) override_args += " --check_spectral_files '${params.msfragger_check_spectral_files}'"
    if (params.msfragger_minimum_peaks != null) override_args += " --minimum_peaks '${params.msfragger_minimum_peaks}'"
    if (params.msfragger_use_topN_peaks != null) override_args += " --use_topN_peaks '${params.msfragger_use_topN_peaks}'"
    if (params.msfragger_min_fragments_modelling != null) override_args += " --min_fragments_modelling '${params.msfragger_min_fragments_modelling}'"
    if (params.msfragger_min_matched_fragments != null) override_args += " --min_matched_fragments '${params.msfragger_min_matched_fragments}'"
    if (params.msfragger_min_sequence_matches != null) override_args += " --min_sequence_matches '${params.msfragger_min_sequence_matches}'"
    if (params.msfragger_minimum_ratio != null) override_args += " --minimum_ratio '${params.msfragger_minimum_ratio}'"
    if (params.msfragger_clear_mz_range != null) override_args += " --clear_mz_range '${params.msfragger_clear_mz_range}'"
    if (params.msfragger_add_Cterm_peptide != null) override_args += " --add_Cterm_peptide '${params.msfragger_add_Cterm_peptide}'"
    if (params.msfragger_add_Nterm_peptide != null) override_args += " --add_Nterm_peptide '${params.msfragger_add_Nterm_peptide}'"
    if (params.msfragger_add_Cterm_protein != null) override_args += " --add_Cterm_protein '${params.msfragger_add_Cterm_protein}'"
    if (params.msfragger_add_Nterm_protein != null) override_args += " --add_Nterm_protein '${params.msfragger_add_Nterm_protein}'"
    if (params.msfragger_add_G_glycine != null) override_args += " --add_G_glycine '${params.msfragger_add_G_glycine}'"
    if (params.msfragger_add_A_alanine != null) override_args += " --add_A_alanine '${params.msfragger_add_A_alanine}'"
    if (params.msfragger_add_S_serine != null) override_args += " --add_S_serine '${params.msfragger_add_S_serine}'"
    if (params.msfragger_add_P_proline != null) override_args += " --add_P_proline '${params.msfragger_add_P_proline}'"
    if (params.msfragger_add_V_valine != null) override_args += " --add_V_valine '${params.msfragger_add_V_valine}'"
    if (params.msfragger_add_T_threonine != null) override_args += " --add_T_threonine '${params.msfragger_add_T_threonine}'"
    if (params.msfragger_add_C_cysteine != null) override_args += " --add_C_cysteine '${params.msfragger_add_C_cysteine}'"
    if (params.msfragger_add_L_leucine != null) override_args += " --add_L_leucine '${params.msfragger_add_L_leucine}'"
    if (params.msfragger_add_I_isoleucine != null) override_args += " --add_I_isoleucine '${params.msfragger_add_I_isoleucine}'"
    if (params.msfragger_add_N_asparagine != null) override_args += " --add_N_asparagine '${params.msfragger_add_N_asparagine}'"
    if (params.msfragger_add_D_aspartic_acid != null) override_args += " --add_D_aspartic_acid '${params.msfragger_add_D_aspartic_acid}'"
    if (params.msfragger_add_Q_glutamine != null) override_args += " --add_Q_glutamine '${params.msfragger_add_Q_glutamine}'"
    if (params.msfragger_add_K_lysine != null) override_args += " --add_K_lysine '${params.msfragger_add_K_lysine}'"
    if (params.msfragger_add_E_glutamic_acid != null) override_args += " --add_E_glutamic_acid '${params.msfragger_add_E_glutamic_acid}'"
    if (params.msfragger_add_M_methionine != null) override_args += " --add_M_methionine '${params.msfragger_add_M_methionine}'"
    if (params.msfragger_add_H_histidine != null) override_args += " --add_H_histidine '${params.msfragger_add_H_histidine}'"
    if (params.msfragger_add_F_phenylalanine != null) override_args += " --add_F_phenylalanine '${params.msfragger_add_F_phenylalanine}'"
    if (params.msfragger_add_R_arginine != null) override_args += " --add_R_arginine '${params.msfragger_add_R_arginine}'"
    if (params.msfragger_add_Y_tyrosine != null) override_args += " --add_Y_tyrosine '${params.msfragger_add_Y_tyrosine}'"
    if (params.msfragger_add_W_tryptophan != null) override_args += " --add_W_tryptophan '${params.msfragger_add_W_tryptophan}'"
    if (params.msfragger_add_B_user_amino_acid != null) override_args += " --add_B_user_amino_acid '${params.msfragger_add_B_user_amino_acid}'"
    if (params.msfragger_add_J_user_amino_acid != null) override_args += " --add_J_user_amino_acid '${params.msfragger_add_J_user_amino_acid}'"
    if (params.msfragger_add_O_user_amino_acid != null) override_args += " --add_O_user_amino_acid '${params.msfragger_add_O_user_amino_acid}'"
    if (params.msfragger_add_U_user_amino_acid != null) override_args += " --add_U_user_amino_acid '${params.msfragger_add_U_user_amino_acid}'"
    if (params.msfragger_add_X_user_amino_acid != null) override_args += " --add_X_user_amino_acid '${params.msfragger_add_X_user_amino_acid}'"
    if (params.msfragger_add_Z_user_amino_acid != null) override_args += " --add_Z_user_amino_acid '${params.msfragger_add_Z_user_amino_acid}'"
    
    """
    generate_msfragger_config.py ${template_file} msfragger_config.params${override_args}
    """
}