// ProSeg crop-level grid search.
// CLI: proseg --xenium <transcripts.parquet> [OPTIONS]
// Transcript assignments are in --output-transcript-metadata (cell_id column).
// Python join step is handled by the downstream PROSEG_ASSIGN_CROP module
// (proseg container has no Python).
process PROSEG_CROP {
    tag "${meta.id}:${crop_id}:${param_hash}"
    label 'process_medium'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    input:
    tuple val(meta), val(crop_id), path(transcripts), val(param_hash),
          val(cell_compactness), val(max_transcript_nucleus_distance),
          val(voxel_size), val(diffusion_probability)

    output:
    tuple val(meta), val(crop_id), val(param_hash),
          path("transcript_metadata.parquet"), path("cell_polygons.geojson.gz"),
          emit: raw

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
    """
    stub:
    """
    touch cell_polygons.geojson.gz transcript_metadata.parquet
    """
}
