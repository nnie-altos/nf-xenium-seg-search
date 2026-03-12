process DETECT_XOA_VERSION {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.0'

    input:
    tuple val(meta), path(experiment_xenium)

    output:
    tuple val(meta), env(XOA_VERSION), emit: version

    script:
    """
    detect_xoa_version.py ${experiment_xenium}
    XOA_VERSION=\$(cat xoa_version.txt)
    """
}
