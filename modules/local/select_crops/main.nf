process SELECT_CROPS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.0'

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
        --n-crops ${n_crops} \\
        --crop-size-um ${crop_size_um} \\
        --pixel-size-um ${pixel_size_um}
    """
}
