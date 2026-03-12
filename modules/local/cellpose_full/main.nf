// Cellpose full-slide segmentation for Stage 2 (optimal parameters).
// Uses lazy tifffile read to avoid OOM on full-slide OME-TIFF.
process CELLPOSE_FULL {
    tag "${meta.id}"
    label 'process_gpu'

    container 'ghcr.io/mouseland/cellpose:3.0.11'

    publishDir "${params.outdir}/${meta.id}/cellpose", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(transcripts), path(morphology_tif),
          val(diameter), val(flow_threshold), val(sharpen_tiff)
    val pixel_size_um

    output:
    tuple val(meta), val("cellpose"),
          path("cells.parquet"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    # Optionally sharpen image before Cellpose (unsharp mask via skimage)
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

    python3 -m cellpose \\
        --image_path \${INPUT_IMG} \\
        --pretrained_model cyto3 \\
        --diameter ${diameter} \\
        --flow_threshold ${flow_threshold} \\
        --no_npy \\
        --save_tif \\
        --savedir .

    python3 - <<'PYEOF'
import numpy as np
import pandas as pd
import tifffile
from pathlib import Path
from skimage.measure import regionprops_table

masks = sorted(Path(".").glob("*_cp_masks.tif"))
if not masks:
    import sys; sys.exit("ERROR: no Cellpose mask file found")
mask_arr = tifffile.imread(str(masks[0]))

props = regionprops_table(mask_arr, properties=["label", "centroid", "area"])
df_cells = pd.DataFrame(props)
df_cells.rename(columns={"label": "cell_id", "centroid-0": "centroid_y_px",
                          "centroid-1": "centroid_x_px"}, inplace=True)
df_cells.to_parquet("cells.parquet", index=False)

# Assign transcripts via pixel-level mask lookup
# Transcripts carry global micron coords; mask pixel 0,0 = slide origin
pixel_size = float("${pixel_size_um}")
t = pd.read_parquet("${transcripts}")
t_px_x = (t["x_location"] / pixel_size).astype(int)
t_px_y = (t["y_location"] / pixel_size).astype(int)
h, w = mask_arr.shape[-2], mask_arr.shape[-1]
valid = (t_px_x >= 0) & (t_px_x < w) & (t_px_y >= 0) & (t_px_y < h)
cell_ids = np.zeros(len(t), dtype=int)
cell_ids[valid.values] = mask_arr[t_px_y[valid].values, t_px_x[valid].values]
t["cell_id"] = cell_ids
t["cell_id"] = t["cell_id"].replace(0, None)
t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
    stub:
    """
    touch cells.parquet transcripts_assigned.parquet
    """
}
