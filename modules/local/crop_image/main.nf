process CROP_IMAGE {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.2'

    input:
    tuple val(meta), path(morphology_tif), path(crops_csv)

    output:
    tuple val(meta), path("*_morphology.tif"), emit: cropped_images

    script:
    """
    python3 - <<'PYEOF'
import csv, subprocess
with open("${crops_csv}") as f:
    for row in csv.DictReader(f):
        subprocess.run([
            "crop_image.py",
            "${morphology_tif}",
            row["px_x_min"], row["px_x_max"],
            row["px_y_min"], row["px_y_max"],
            row["crop_id"],
        ], check=True)
PYEOF
    """
    stub:
    """
    touch stub_crop1_morphology.tif
    """
}
