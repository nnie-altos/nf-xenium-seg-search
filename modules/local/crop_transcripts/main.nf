process CROP_TRANSCRIPTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.2'

    input:
    tuple val(meta), path(transcripts), path(crops_csv)

    output:
    tuple val(meta), path("*_transcripts.parquet"), emit: cropped_transcripts

    script:
    """
    python3 - <<'PYEOF'
import csv, subprocess
with open("${crops_csv}") as f:
    for row in csv.DictReader(f):
        subprocess.run([
            "crop_transcripts.py",
            "${transcripts}",
            row["x_min_um"], row["x_max_um"],
            row["y_min_um"], row["y_max_um"],
            row["crop_id"],
        ], check=True)
PYEOF
    """
    stub:
    """
    touch stub_crop1_transcripts.parquet
    """
}
