// SEGGER full-slide segmentation for Stage 2 (optimal parameters).
process SEGGER_FULL {
    tag "${meta.id}"
    label 'process_gpu'

    container 'docker.io/altoslabscom/xenium-processing-gpu:segger-1.0.14'

    publishDir "${params.outdir}/${meta.id}/segger", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(xenium_bundle, stageAs: "bundle/"),
          val(tile_size), val(min_transcripts_per_cell)
    path segger_model

    output:
    tuple val(meta), val("segger"),
          path("segger_cells.parquet"),
          path("transcripts_assigned.parquet"),
          emit: results

    script:
    def CREATE_SCRIPT  = "/workspace/segger_dev/src/segger/cli/create_dataset_fast.py"
    def PREDICT_SCRIPT = "/workspace/segger_dev/src/segger/cli/predict_fast.py"
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    # ── Prepare bundle ────────────────────────────────────────────────────────────
    python3 - <<'PYEOF'
import pyarrow.parquet as pq, os

os.makedirs('bundle_stats', exist_ok=True)
for f in ['transcripts.parquet', 'nucleus_boundaries.parquet']:
    src = f'bundle/{f}'
    if not os.path.exists(src):
        continue
    t = pq.read_table(src)
    pq.write_table(t, f'bundle_stats/{f}', write_statistics=True, compression='snappy')

for item in os.listdir('bundle'):
    dst = f'bundle_stats/{item}'
    if not os.path.exists(dst):
        os.symlink(os.path.realpath(f'bundle/{item}'), dst)

if not os.path.exists('bundle_stats/nucleus_boundaries.parquet') and \
        os.path.exists('bundle_stats/cell_boundaries.parquet'):
    os.symlink(os.path.realpath('bundle_stats/cell_boundaries.parquet'),
               'bundle_stats/nucleus_boundaries.parquet')
PYEOF

    # ── Create SEGGER dataset ─────────────────────────────────────────────────────
    python3 ${CREATE_SCRIPT} \\
        --base_dir    bundle_stats \\
        --data_dir    segger_dataset \\
        --sample_type xenium \\
        --tile_width  ${tile_size} \\
        --tile_height ${tile_size} \\
        --n_workers   ${task.cpus}

    # ── Run SEGGER prediction ─────────────────────────────────────────────────────
    GPU_IDS=\$(python3 -c "
import torch
n = torch.cuda.device_count()
print(','.join(str(i) for i in range(n)) if n > 0 else '0')
" 2>/dev/null || echo "0")

    python3 ${PREDICT_SCRIPT} \\
        --models_dir       ${segger_model} \\
        --segger_data_dir  segger_dataset \\
        --transcripts_file bundle/transcripts.parquet \\
        --benchmarks_dir   benchmarks_dir \\
        --batch_size       4 \\
        --use_cc           False \\
        --knn_method       kd_tree \\
        --num_workers      ${task.cpus} \\
        --gpu_ids          \${GPU_IDS}

    # ── Post-process ──────────────────────────────────────────────────────────────
    python3 - <<'PYEOF'
import pandas as pd, sys
from pathlib import Path

tx_files = sorted(Path("benchmarks_dir").glob("*/segger_transcripts.parquet"))
if not tx_files:
    sys.exit("ERROR: SEGGER produced no transcript output in benchmarks_dir")

t = pd.read_parquet(str(tx_files[0]))

if "segger_cell_id" in t.columns:
    t = t.rename(columns={"segger_cell_id": "cell_id"})
elif "cell_id" not in t.columns:
    cell_col = [c for c in t.columns if "cell" in c.lower()]
    if cell_col:
        t = t.rename(columns={cell_col[0]: "cell_id"})
    else:
        t["cell_id"] = None

min_tx = int("${min_transcripts_per_cell}")
if min_tx > 0 and t["cell_id"].notna().any():
    counts = t.dropna(subset=["cell_id"]).groupby("cell_id").size()
    keep_cells = set(counts[counts >= min_tx].index)
    t["cell_id"] = t["cell_id"].where(t["cell_id"].isin(keep_cells), other=None)

t.to_parquet("transcripts_assigned.parquet", index=False)

assigned = t.dropna(subset=["cell_id"])
cells = assigned.groupby("cell_id").size().reset_index(name="transcript_count")
cells.to_parquet("segger_cells.parquet", index=False)
PYEOF
    """
    stub:
    """
    touch segger_cells.parquet transcripts_assigned.parquet
    """
}
