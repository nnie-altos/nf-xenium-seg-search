// Cellpose full-slide segmentation — produces mask TIF for downstream Baysor prior.
// Stage 2 only; uses optimal Cellpose parameters from Stage 1.
process CELLPOSE_MASK_FULL {
    tag "${meta.id}"
    label 'process_gpu'

    container 'ghcr.io/mouseland/cellpose:3.0.11'

    input:
    tuple val(meta), path(morphology_tif),
          val(diameter), val(flow_threshold), val(sharpen_tiff)
    val nucleus_segmentation_only

    output:
    tuple val(meta), path("cp_mask.tif"), emit: mask

    script:
    """
    if [ "${sharpen_tiff}" = "true" ]; then
        python3 - <<'SHARP'
import tifffile, numpy as np
from skimage.filters import unsharp_mask
img = tifffile.imread("${morphology_tif}")
sharp = unsharp_mask(img.astype(np.float64), radius=2.0, amount=1.5)
if img.dtype == np.uint16:
    tifffile.imwrite("input_for_cellpose.tif",
                     np.clip(sharp * 65535, 0, 65535).astype(np.uint16))
else:
    tifffile.imwrite("input_for_cellpose.tif",
                     np.clip(sharp * 255, 0, 255).astype(np.uint8))
SHARP
        INPUT_IMG="input_for_cellpose.tif"
    else
        INPUT_IMG="${morphology_tif}"
    fi

    CELLPOSE_MODEL="${nucleus_segmentation_only}" == "true" ? "nuclei" : "cyto3"
    python3 -m cellpose \\
        --image_path \${INPUT_IMG} \\
        --pretrained_model \${CELLPOSE_MODEL} \\
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
