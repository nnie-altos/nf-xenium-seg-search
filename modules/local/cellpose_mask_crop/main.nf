// Cellpose crop segmentation — produces mask TIF for downstream Baysor prior.
// Runs once per crop; Baysor param combos reuse this mask.
process CELLPOSE_MASK_CROP {
    tag "${meta.id}:${crop_id}"
    label 'process_gpu'

    container 'docker.io/altoslabscom/cellpose:4.0.8'

    input:
    tuple val(meta), val(crop_id), path(image)
    val nucleus_segmentation_only
    val diameter
    val flow_threshold

    output:
    tuple val(meta), val(crop_id), path("cp_mask.tif"), emit: mask

    script:
    def cellpose_model = nucleus_segmentation_only.toString() == "true" ? "nuclei" : "cyto3"
    """
    python3 -m cellpose \\
        --image_path ${image} \\
        --pretrained_model ${cellpose_model} \\
        --diameter ${diameter} \\
        --flow_threshold ${flow_threshold} \\
        --no_npy \\
        --save_tif \\
        --savedir .

    mv *_cp_masks.tif cp_mask.tif
    """
    stub:
    """
    touch cp_mask.tif
    """
}
