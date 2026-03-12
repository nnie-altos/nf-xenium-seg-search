// SEGGER full-slide segmentation for Stage 2 (optimal parameters).
process SEGGER_FULL {
    tag "${meta.id}"
    label 'process_gpu'

    container 'danielunyi42/segger_cuda118_py311:latest'

    publishDir "${params.outdir}/${meta.id}/segger", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(transcripts),
          val(min_transcripts_per_cell), val(tile_size)
    path segger_model

    output:
    tuple val(meta), val("segger"),
          path("segger_cells.parquet"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    """
    segger predict \\
        --transcript-file ${transcripts} \\
        --model-dir       ${segger_model} \\
        --output-dir      . \\
        --min-transcripts ${min_transcripts_per_cell} \\
        --tile-size       ${tile_size}

    CELLS_OUT=\$(ls *_cells.parquet 2>/dev/null | head -1)
    TRANS_OUT=\$(ls *_transcripts.parquet 2>/dev/null | head -1)

    [ -n "\$CELLS_OUT" ] && cp "\$CELLS_OUT" segger_cells.parquet || \\
        python3 -c "import pandas as pd; pd.DataFrame(columns=['cell_id']).to_parquet('segger_cells.parquet', index=False)"

    [ -n "\$TRANS_OUT" ] && cp "\$TRANS_OUT" transcripts_assigned.parquet || \\
        python3 -c "
import pandas as pd
t = pd.read_parquet('${transcripts}')
if 'cell_id' not in t.columns: t['cell_id'] = None
t.to_parquet('transcripts_assigned.parquet', index=False)
"
    """
}
