process CRYSTAL_C_CONFIG {
    publishDir "${output_dir}/02-Crystal-C/",
        mode: params.publish_dir_mode,
        pattern: "crystalc_config.params"

    input:
    val(output_dir)
    path(template_file)
    path(database_fasta)

    output:
    path("crystalc_config.params"), emit: crystal_c_config

    script:
    def database_path = database_fasta.toString()
    
    // Build override arguments
    def override_args = " --thread ${task.cpus} --fasta '${database_path}' --output_location '.'"
    if (params.crystal_c_precursor_charge != null) override_args += " --precursor_charge '${params.crystal_c_precursor_charge}'"
    if (params.crystal_c_isotope_number != null) override_args += " --isotope_number '${params.crystal_c_isotope_number}'"
    if (params.crystal_c_precursor_mass != null) override_args += " --precursor_mass '${params.crystal_c_precursor_mass}'"
    if (params.crystal_c_precursor_isolation_window != null) override_args += " --precursor_isolation_window '${params.crystal_c_precursor_isolation_window}'"
    if (params.crystal_c_correct_isotope_error != null) override_args += " --correct_isotope_error '${params.crystal_c_correct_isotope_error}'"
    if (params.crystal_c_raw_file_location != null) override_args += " --raw_file_location '${params.crystal_c_raw_file_location}'"
    if (params.crystal_c_raw_file_extension != null) override_args += " --raw_file_extension '${params.crystal_c_raw_file_extension}'"
    
    """
    generate_crystalc_config.py ${template_file} crystalc_config.params${override_args}
    """
}

process CRYSTAL_C {
    publishDir "${output_dir}/02-Crystal-C/",
        mode: params.publish_dir_mode,
        pattern: "*.pepXML"
    
    input:
    val(output_dir)
    path(pepxml_file)
    path(crystal_c_config)

    output:
    path("*.pepXML"), emit: pepxml_recalibrated
    path("versions.yml"), emit: versions

    script:
    def memory_gb = task.memory.toGiga()
    def config_file = crystal_c_config.toString()
    """
    # Crystal-C mass recalibration for open searches
    java -Xmx${memory_gb}G -cp "CrystalC-1.2.1.jar" crystalc.Run \\
        ${config_file} \\
        ${pepxml_file}
    
    # Version info
    echo '"${task.process}":' > versions.yml
    echo "    crystal-c: \$(java -cp CrystalC-1.2.1.jar crystalc.Run 2>&1 | head -n1 || echo 'unknown')" >> versions.yml
    echo "    java: \$(java -version 2>&1 | head -n1 | cut -d' ' -f3 | tr -d '\"')" >> versions.yml
    """
}

