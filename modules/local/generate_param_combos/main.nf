process GENERATE_PARAM_COMBOS {
    tag "param_combos"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.0'

    input:
    path param_grids_yaml

    output:
    path "param_manifest.json", emit: manifest
    path "*_params.tsv",        emit: param_tsvs

    script:
    """
    generate_param_combos.py ${param_grids_yaml}
    """
    stub:
    """
    echo '{}' > param_manifest.json
    printf 'param_hash\tcell_compactness\tmax_transcript_nucleus_distance\tvoxel_size\tdiffusion_probability\nstub_p1\t0.04\t60\t1.0\t0.2\n' > proseg_params.tsv
    printf 'param_hash\tdiameter\tflow_threshold\tsharpen_tiff\nstub_c1\t30\t0.4\tfalse\n' > cellpose_params.tsv
    printf 'param_hash\tmin_transcripts_per_cell\ttile_size\nstub_s1\t5\t120\n' > segger_params.tsv
    printf 'param_hash\texpansion_distance\tdapi_filter\tboundary_stain\nstub_x1\t5\t100\ttrue\n' > xr_params.tsv
    """
}
