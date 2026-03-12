// Cellpose crop-level grid search.
// Runs Cellpose on the cropped OME-TIFF, then assigns transcripts to cells
// by looking up each transcript's pixel coordinate in the Cellpose mask array.
// crops_csv is used to read the crop bounding box (x_min_um, y_min_um) for
// converting micron transcript coords → pixel coords within the crop image.
process CELLPOSE_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_gpu'

    container 'ghcr.io/mouseland/cellpose:3.0.11'

    input:
    tuple val(meta), val(crop_id), path(transcripts), path(image), path(crops_csv),
          val(param_hash), val(diameter), val(flow_threshold), val(sharpen_tiff)
    val pixel_size_um

    output:
    tuple val(meta), val(crop_id), val("cellpose"), val(param_hash),
          path("cells.parquet"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    # Optionally sharpen image before Cellpose (unsharp mask via skimage)
    if [ "${sharpen_tiff}" = "true" ]; then
        python3 - <<'SHARP'
import tifffile, numpy as np
from skimage.filters import unsharp_mask
img = tifffile.imread("${image}")
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
        INPUT_IMG="${image}"
    fi

    # Run Cellpose segmentation
    python3 -m cellpose \\
        --image_path \${INPUT_IMG} \\
        --pretrained_model cyto3 \\
        --diameter ${diameter} \\
        --flow_threshold ${flow_threshold} \\
        --no_npy \\
        --save_tif \\
        --savedir .

    # Post-process: build cells.parquet and assign transcripts via mask lookup
    python3 - <<'PYEOF'
import numpy as np
import pandas as pd
import tifffile
from pathlib import Path
from skimage.measure import regionprops_table

# ── Find Cellpose mask output ─────────────────────────────────────────────────
masks = sorted(Path(".").glob("*_cp_masks.tif"))
if not masks:
    import sys; sys.exit("ERROR: no Cellpose mask file found")
mask_arr = tifffile.imread(str(masks[0]))  # (Y, X) int array; 0 = background

# ── Build cells.parquet from mask properties ──────────────────────────────────
props = regionprops_table(mask_arr, properties=["label", "centroid", "area"])
df_cells = pd.DataFrame(props)
df_cells.rename(columns={"label": "cell_id", "centroid-0": "centroid_y_px",
                          "centroid-1": "centroid_x_px"}, inplace=True)
df_cells.to_parquet("cells.parquet", index=False)

# ── Assign transcripts via mask array ─────────────────────────────────────────
crops = pd.read_csv("${crops_csv}")
row = crops[crops["crop_id"] == "${crop_id}"].iloc[0]
x_min_um = float(row["x_min_um"])
y_min_um = float(row["y_min_um"])
pixel_size = float("${pixel_size_um}")

t = pd.read_parquet("${transcripts}")
# Convert micron coords to pixel coords within the crop image
t_px_x = ((t["x_location"] - x_min_um) / pixel_size).astype(int)
t_px_y = ((t["y_location"] - y_min_um) / pixel_size).astype(int)

h, w = mask_arr.shape[-2], mask_arr.shape[-1]
valid = (t_px_x >= 0) & (t_px_x < w) & (t_px_y >= 0) & (t_px_y < h)
cell_ids = np.zeros(len(t), dtype=int)
cell_ids[valid.values] = mask_arr[t_px_y[valid].values, t_px_x[valid].values]

t["cell_id"] = cell_ids
t["cell_id"] = t["cell_id"].replace(0, None)  # 0 = background → unassigned
t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
    stub:
    """
    touch cells.parquet transcripts_assigned.parquet
    """
}
