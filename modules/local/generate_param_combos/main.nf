process GENERATE_PARAM_COMBOS {
    tag "param_combos"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.1'

    input:
    path param_grids_yaml
    val  search_strategy
    val  nucleus_segmentation_only

    output:
    path "param_manifest.json", emit: manifest
    path "*_params.tsv",        emit: param_tsvs

    script:
    def nucleus_flag = nucleus_segmentation_only.toString() == 'true' ? '--nucleus-segmentation-only' : ''
    """
    generate_param_combos.py ${param_grids_yaml} --strategy ${search_strategy} ${nucleus_flag}
    """
    stub:
    """
    echo '{}' > param_manifest.json
    printf 'param_hash\tcell_compactness\tmax_transcript_nucleus_distance\tvoxel_size\tdiffusion_probability\nstub_p1\t0.04\t60\t1.0\t0.2\n' > proseg_params.tsv
    printf 'param_hash\tdiameter\tflow_threshold\tsharpen_tiff\nstub_c1\t30\t0.4\tfalse\n' > cellpose_params.tsv
    printf 'param_hash\tmin_transcripts_per_cell\td_max\ttile_size\nstub_s1\t25\t2.0\t120\n' > segger_params.tsv
    printf 'param_hash\tprior_segmentation_confidence\tmin_molecules_per_cell\nstub_b1\t0.5\t50\n' > cellpose_baysor_params.tsv
    """
}
