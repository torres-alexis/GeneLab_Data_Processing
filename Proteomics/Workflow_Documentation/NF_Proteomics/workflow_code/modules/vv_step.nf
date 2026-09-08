process VV_STEP {
    tag "${name}"

    input:
        val name
        path files

    output:
        path("vv_${name}.ok"), emit: ok

    script:
    """
    vv_step.py --name ${name} ${files}
    touch vv_${name}.ok
    """
}
