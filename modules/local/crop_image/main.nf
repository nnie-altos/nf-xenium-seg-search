process CROP_IMAGE {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/altos-labs/xenium-seg-search:0.1.0'

    input:
    tuple val(meta), path(morphology_tif), path(crops_csv)

    output:
    tuple val(meta), path("*_morphology.tif"), emit: cropped_images

    script:
    """
    crop_image.py ${crops_csv} ${morphology_tif}
    """
}
