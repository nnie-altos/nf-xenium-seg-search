// ProSeg full-slide segmentation for Stage 2 (optimal parameters).
process PROSEG_FULL {
    tag "${meta.id}"
    label 'process_high'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    publishDir "${params.outdir}/${meta.id}/proseg", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(transcripts),
          val(cell_perimeter_ratio_threshold), val(cell_aspect_ratio_limit),
          val(cell_size_min), val(cell_size_max), val(nuclei_distance_threshold)

    output:
    tuple val(meta), val("proseg"),
          path("cell_polygons.geojson.gz"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    proseg xenium ${transcripts} \\
        --cell-polygons              cell_polygons.geojson.gz \\
        --transcript-assignments     transcript_assignments.csv \\
        --cell-metadata              cell_metadata.csv \\
        --expected-counts            expected_counts.csv \\
        --cell-perimeter-ratio-threshold ${cell_perimeter_ratio_threshold} \\
        --cell-aspect-ratio-limit        ${cell_aspect_ratio_limit} \\
        --cell-size-min                  ${cell_size_min} \\
        --cell-size-max                  ${cell_size_max} \\
        --nuclei-distance-threshold      ${nuclei_distance_threshold}

    python3 - <<'PYEOF'
import pandas as pd, sys
t = pd.read_parquet("${transcripts}")
a = pd.read_csv("transcript_assignments.csv")
if "transcript_id" in a.columns and "transcript_id" in t.columns:
    t = t.merge(a[["transcript_id", "cell_id"]], on="transcript_id", how="left")
elif len(a) == len(t):
    t["cell_id"] = a.iloc[:, -1].values
else:
    print("WARNING: transcript assignment join failed", file=sys.stderr)
    t["cell_id"] = None
t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
}
