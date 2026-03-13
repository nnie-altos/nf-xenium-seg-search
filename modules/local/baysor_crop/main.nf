// Baysor crop-level grid search using Cellpose mask as prior.
// Only prior_segmentation_confidence and min_molecules_per_cell are varied;
// Cellpose mask is shared across all Baysor param combos for the same crop.
//
// Transcript coordinates are converted to crop-local pixels (matching the mask TIF
// coordinate space) before running Baysor. Results are mapped back to the original
// parquet by index.
process BAYSOR_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_high'

    container 'khersameesh24/baysor:0.7.1'

    input:
    tuple val(meta), val(crop_id), path(transcripts), path(mask),
          path(crops_csv), val(param_hash),
          val(prior_segmentation_confidence), val(min_molecules_per_cell)
    val pixel_size_um
    val scale   // Cellpose diameter in pixels (= Baysor scale)

    output:
    tuple val(meta), val(crop_id), val("cellpose_baysor"), val(param_hash),
          path("cells.parquet"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    # ── Write Baysor TOML (pixel-space coordinates, no z) ────────────────────
    cat > baysor.toml << 'TOML'
[data]
x = "x_px"
y = "y_px"
gene = "feature_name"
min_molecules_per_gene = 1
min_molecules_per_cell = ${min_molecules_per_cell}

[segmentation]
unassigned_prior_label = "UNASSIGNED"
prior_segmentation_confidence = ${prior_segmentation_confidence}

[plotting]
min_pixels_per_cell = 5
TOML

    # ── Convert transcripts to crop-local pixel coordinates ──────────────────
    python3 - <<'PYEOF'
import pandas as pd

crops = pd.read_csv("${crops_csv}")
row = crops[crops["crop_id"] == "${crop_id}"].iloc[0]
x_min_um   = float(row["x_min_um"])
y_min_um   = float(row["y_min_um"])
pixel_size = float("${pixel_size_um}")

t = pd.read_parquet("${transcripts}").reset_index(drop=True)
t["x_px"] = (t["x_location"] - x_min_um) / pixel_size
t["y_px"] = (t["y_location"] - y_min_um) / pixel_size
t.to_csv("transcripts_px.csv", index=True, index_label="__orig_idx")
PYEOF

    # ── Run Baysor ────────────────────────────────────────────────────────────
    baysor run \\
        transcripts_px.csv \\
        ${mask} \\
        --scale=${scale} \\
        --output=segmentation.csv \\
        --config=baysor.toml \\
        --polygon-format=GeometryCollectionLegacy

    # ── Post-process: standard cells.parquet + transcripts_assigned.parquet ──
    python3 - <<'PYEOF'
import pandas as pd, sys

seg = pd.read_csv("segmentation.csv")
if "cell" not in seg.columns:
    sys.exit("ERROR: Baysor segmentation.csv missing 'cell' column")

# Map Baysor cell assignments back to original transcript rows by index
t_orig = pd.read_parquet("${transcripts}").reset_index(drop=True)
idx_map = seg.set_index("__orig_idx")["cell"]
t_orig["cell_id"] = t_orig.index.map(idx_map)
# Baysor uses "0" (string) for unassigned in some versions, empty string in others
t_orig["cell_id"] = t_orig["cell_id"].replace({"0": None, "": None, 0: None})
t_orig.to_parquet("transcripts_assigned.parquet", index=False)

assigned = t_orig.dropna(subset=["cell_id"])
if assigned.empty:
    import pandas as pd
    cells = pd.DataFrame(columns=["cell_id", "centroid_x_um", "centroid_y_um"])
else:
    cells = (
        assigned
        .groupby("cell_id")[["x_location", "y_location"]]
        .mean()
        .rename(columns={"x_location": "centroid_x_um", "y_location": "centroid_y_um"})
        .reset_index()
    )
cells.to_parquet("cells.parquet", index=False)
print(f"Baysor: {len(cells)} cells, {len(assigned)} assigned transcripts",
      file=__import__('sys').stderr)
PYEOF
    """
    stub:
    """
    touch cells.parquet transcripts_assigned.parquet
    """
}
