// ProSeg crop-level grid search process.
// ProSeg CLI: `proseg xenium <transcripts.parquet> [OPTIONS]`
// Outputs cell_polygons.geojson.gz and transcript_assignments.csv;
// the post-step joins assignments back to transcripts for scoring.
process PROSEG_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_medium'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    input:
    tuple val(meta), val(crop_id), path(transcripts), val(param_hash),
          val(cell_perimeter_ratio_threshold), val(cell_aspect_ratio_limit),
          val(cell_size_min), val(cell_size_max), val(nuclei_distance_threshold)

    output:
    tuple val(meta), val(crop_id), val("proseg"), val(param_hash),
          path("cell_polygons.geojson.gz"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    proseg xenium ${transcripts} \\
        --cell-polygons         cell_polygons.geojson.gz \\
        --transcript-assignments transcript_assignments.csv \\
        --cell-metadata         cell_metadata.csv \\
        --cell-perimeter-ratio-threshold ${cell_perimeter_ratio_threshold} \\
        --cell-aspect-ratio-limit        ${cell_aspect_ratio_limit} \\
        --cell-size-min                  ${cell_size_min} \\
        --cell-size-max                  ${cell_size_max} \\
        --nuclei-distance-threshold      ${nuclei_distance_threshold}

    # Join ProSeg transcript assignments back to input parquet
    python3 - <<'PYEOF'
import pandas as pd, sys
t = pd.read_parquet("${transcripts}")
a = pd.read_csv("transcript_assignments.csv")
# ProSeg outputs: transcript_id (index), cell_id
if "transcript_id" in a.columns and "transcript_id" in t.columns:
    t = t.merge(a[["transcript_id", "cell_id"]], on="transcript_id", how="left")
elif len(a) == len(t):
    t["cell_id"] = a.iloc[:, -1].values
else:
    print("WARNING: transcript assignment join failed; marking all unassigned", file=sys.stderr)
    t["cell_id"] = None
t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
}
