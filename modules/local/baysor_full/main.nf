// Baysor full-slide segmentation using Cellpose mask as prior.
// Stage 2 only; uses optimal cellpose_baysor parameters from Stage 1.
//
// Transcript coordinates are converted to pixel space to match the Cellpose mask TIF.
// Results are mapped back to the original parquet by index.
process BAYSOR_FULL {
    tag "${meta.id}"
    label 'process_high'

    container 'khersameesh24/baysor:0.7.1'

    publishDir "${params.outdir}/${meta.id}/cellpose_baysor", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(transcripts), path(mask),
          val(prior_segmentation_confidence), val(min_molecules_per_cell),
          val(scale)   // optimal Cellpose diameter in pixels
    val pixel_size_um

    output:
    tuple val(meta), val("cellpose_baysor"),
          path("cells.parquet"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    # ── Write Baysor TOML (pixel-space coordinates) ───────────────────────────
    cat > baysor.toml << 'TOML'
[data]
x = "x_px"
y = "y_px"
gene = "feature_name"
min_molecules_per_gene = 10
exclude_genes = "NegControl*,BLANK_*,antisense_*"
min_molecules_per_cell = ${min_molecules_per_cell}

[segmentation]
unassigned_prior_label = "UNASSIGNED"
prior_segmentation_confidence = ${prior_segmentation_confidence}

[plotting]
min_pixels_per_cell = 10
TOML

    # ── Convert transcripts to pixel coordinates ──────────────────────────────
    python3 - <<'PYEOF'
import pandas as pd

pixel_size = float("${pixel_size_um}")
t = pd.read_parquet("${transcripts}").reset_index(drop=True)
t["x_px"] = t["x_location"] / pixel_size
t["y_px"] = t["y_location"] / pixel_size
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

t_orig = pd.read_parquet("${transcripts}").reset_index(drop=True)
idx_map = seg.set_index("__orig_idx")["cell"]
t_orig["cell_id"] = t_orig.index.map(idx_map)
t_orig["cell_id"] = t_orig["cell_id"].replace({"0": None, "": None, 0: None})
t_orig.to_parquet("transcripts_assigned.parquet", index=False)

assigned = t_orig.dropna(subset=["cell_id"])
if assigned.empty:
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
