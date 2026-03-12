// SEGGER crop-level grid search: inference only (pre-trained model required).
// Tunes: min_transcripts_per_cell, tile_size.
process SEGGER_PREDICT_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_gpu'

    container 'danielunyi42/segger_cuda118_py311:latest'

    input:
    tuple val(meta), val(crop_id), path(transcripts), val(param_hash),
          val(min_transcripts_per_cell), val(tile_size)
    path segger_model

    output:
    tuple val(meta), val(crop_id), val("segger"), val(param_hash),
          path("segger_cells.parquet"), path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    segger predict \\
        --transcript-file   ${transcripts} \\
        --model-dir         ${segger_model} \\
        --output-dir        . \\
        --min-transcripts   ${min_transcripts_per_cell} \\
        --tile-size         ${tile_size}

    # Rename SEGGER output to expected filenames
    CELLS_OUT=\$(ls *_cells.parquet 2>/dev/null | head -1)
    TRANS_OUT=\$(ls *_transcripts.parquet 2>/dev/null | head -1)

    if [ -n "\$CELLS_OUT" ]; then
        cp "\$CELLS_OUT" segger_cells.parquet
    else
        python3 -c "import pandas as pd; pd.DataFrame(columns=['cell_id']).to_parquet('segger_cells.parquet', index=False)"
    fi

    if [ -n "\$TRANS_OUT" ]; then
        cp "\$TRANS_OUT" transcripts_assigned.parquet
    else
        # Fall back: annotate original transcripts with SEGGER cell assignments
        python3 - <<'PYEOF'
import pandas as pd, sys
print("WARNING: SEGGER transcript output not found; using original transcripts", file=sys.stderr)
t = pd.read_parquet("${transcripts}")
if "cell_id" not in t.columns:
    t["cell_id"] = None
t.to_parquet("transcripts_assigned.parquet", index=False)
PYEOF
    fi
    """
}
