// Joins ProSeg transcript metadata (cell_id) back to the original transcripts parquet.
// Runs in the xenium-seg-search container because the proseg container has no Python.
process PROSEG_ASSIGN_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.4'

    input:
    tuple val(meta), val(crop_id), val(param_hash),
          path(transcripts), path(transcript_metadata), path(cell_polygons)

    output:
    tuple val(meta), val(crop_id), val("proseg"), val(param_hash),
          path("cell_polygons.geojson.gz"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    python3 - <<'PYEOF'
import pandas as pd, sys

t = pd.read_parquet("${transcripts}")
m = pd.read_parquet("${transcript_metadata}")

print(f"  transcript_metadata columns: {list(m.columns)}", file=sys.stderr)

# ProSeg uses 'assignment' not 'cell_id'; normalise to cell_id
assign_col = None
for cname in ["cell_id", "assignment"]:
    if cname in m.columns:
        assign_col = cname
        break
if assign_col is None:
    # Last resort: any column with 'cell' or 'assign' in the name
    for c in m.columns:
        if "cell" in c.lower() or "assign" in c.lower():
            assign_col = c
            break

if assign_col is None:
    print("WARNING: could not find cell assignment column in transcript metadata", file=sys.stderr)
    t["cell_id"] = None
elif "transcript_id" in m.columns and "transcript_id" in t.columns:
    merged = m[["transcript_id", assign_col]].rename(columns={assign_col: "cell_id"})
    t = t.merge(merged, on="transcript_id", how="left")
elif len(m) == len(t):
    t = t.copy()
    t["cell_id"] = m[assign_col].values
else:
    print("WARNING: transcript_metadata length mismatch; skipping assignment", file=sys.stderr)
    t["cell_id"] = None

# ProSeg marks unassigned transcripts with 0 or empty string
if "cell_id" in t.columns:
    t["cell_id"] = t["cell_id"].replace(0, None).replace("", None)

t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    """
    stub:
    """
    touch cell_polygons.geojson.gz transcripts_assigned.parquet
    """
}
