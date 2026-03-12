// ProSeg crop-level grid search.
// CLI: proseg --xenium <transcripts.parquet> [OPTIONS]
// Transcript assignments are in --output-transcript-metadata (cell_id column).
process PROSEG_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_medium'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    input:
    tuple val(meta), val(crop_id), path(transcripts), val(param_hash),
          val(cell_compactness), val(max_transcript_nucleus_distance),
          val(voxel_size), val(diffusion_probability)

    output:
    tuple val(meta), val(crop_id), val("proseg"), val(param_hash),
          path("cell_polygons.geojson.gz"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    proseg --xenium ${transcripts} \\
        --output-cell-polygons      cell_polygons.geojson.gz \\
        --output-transcript-metadata transcript_metadata.parquet \\
        --output-cell-metadata      cell_metadata.parquet \\
        --cell-compactness          ${cell_compactness} \\
        --max-transcript-nucleus-distance ${max_transcript_nucleus_distance} \\
        --voxel-size                ${voxel_size} \\
        --diffusion-probability     ${diffusion_probability} \\
        --ignore-z-coord

    # Join ProSeg transcript metadata (cell_id) back to original transcripts parquet
    python3 - <<'PYEOF'
import pandas as pd, sys

t = pd.read_parquet("${transcripts}")
m = pd.read_parquet("transcript_metadata.parquet")

print(f"  transcript_metadata columns: {list(m.columns)}", file=sys.stderr)

# ProSeg transcript metadata contains a cell_id column; join by index or transcript_id
if "transcript_id" in m.columns and "transcript_id" in t.columns:
    t = t.merge(m[["transcript_id", "cell_id"]], on="transcript_id", how="left")
elif "cell_id" in m.columns and len(m) == len(t):
    t = t.copy()
    t["cell_id"] = m["cell_id"].values
else:
    # Try positional — last resort
    cell_col = [c for c in m.columns if "cell" in c.lower()]
    if cell_col:
        t = t.copy()
        t["cell_id"] = m[cell_col[0]].values
    else:
        print("WARNING: could not find cell assignment column in transcript metadata", file=sys.stderr)
        t["cell_id"] = None

t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
}
