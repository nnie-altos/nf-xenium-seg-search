// ProSeg full-slide segmentation for Stage 2 (optimal parameters).
// Python join step is handled by the downstream PROSEG_ASSIGN_FULL module
// (proseg container has no Python).
process PROSEG_FULL {
    tag "${meta.id}"
    label 'process_high'

    container 'ghcr.io/dcjones/proseg:v3.1.0'

    input:
    tuple val(meta), path(transcripts),
          val(cell_compactness), val(max_transcript_nucleus_distance),
          val(voxel_size), val(diffusion_probability)

    output:
    tuple val(meta), path("transcript_metadata.parquet"), path("cell_polygons.geojson.gz"),
          emit: raw

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
    """
    stub:
    """
    touch cell_polygons.geojson.gz transcript_metadata.parquet
    """
}
