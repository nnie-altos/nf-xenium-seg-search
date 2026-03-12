// ProSeg full-slide segmentation for Stage 2 (optimal parameters).
process PROSEG_FULL {
    tag "${meta.id}"
    label 'process_high'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    publishDir "${params.outdir}/${meta.id}/proseg", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(transcripts),
          val(cell_compactness), val(max_transcript_nucleus_distance),
          val(voxel_size), val(diffusion_probability)

    output:
    tuple val(meta), val("proseg"),
          path("cell_polygons.geojson.gz"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    proseg --xenium ${transcripts} \\
        --output-cell-polygons      cell_polygons.geojson.gz \\
        --output-transcript-metadata transcript_metadata.parquet \\
        --output-cell-metadata      cell_metadata.parquet \\
        --output-expected-counts    expected_counts.parquet \\
        --cell-compactness          ${cell_compactness} \\
        --max-transcript-nucleus-distance ${max_transcript_nucleus_distance} \\
        --voxel-size                ${voxel_size} \\
        --diffusion-probability     ${diffusion_probability} \\
        --ignore-z-coord

    python3 - <<'PYEOF'
import pandas as pd, sys

t = pd.read_parquet("${transcripts}")
m = pd.read_parquet("transcript_metadata.parquet")

print(f"  transcript_metadata columns: {list(m.columns)}", file=sys.stderr)

if "transcript_id" in m.columns and "transcript_id" in t.columns:
    t = t.merge(m[["transcript_id", "cell_id"]], on="transcript_id", how="left")
elif "cell_id" in m.columns and len(m) == len(t):
    t = t.copy()
    t["cell_id"] = m["cell_id"].values
else:
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
