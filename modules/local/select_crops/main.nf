process SELECT_CROPS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.1'

    publishDir "${params.outdir}/${meta.id}/crops", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(nucleus_boundaries)
    val n_crops
    val crop_size_um
    val pixel_size_um

    output:
    tuple val(meta), path("crops.csv"), emit: crops

    script:
    """
    select_crops.py \\
        ${nucleus_boundaries} \\
        ${n_crops} \\
        ${crop_size_um} \\
        ${pixel_size_um}
    """
    stub:
    """
    printf 'crop_id,x_min_um,y_min_um,x_max_um,y_max_um\nstub_crop1,0,0,500,500\n' > crops.csv
    """
}
