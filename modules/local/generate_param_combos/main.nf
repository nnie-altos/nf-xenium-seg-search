process GENERATE_PARAM_COMBOS {
    tag "param_combos"
    label 'process_low'

    container 'ghcr.io/altos-labs/xenium-seg-search:0.1.0'

    input:
    path param_grids_yaml

    output:
    path "param_manifest.json", emit: manifest
    path "*_params.tsv",        emit: param_tsvs

    script:
    """
    generate_param_combos.py ${param_grids_yaml}
    """
}
