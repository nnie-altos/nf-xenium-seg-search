// Crop transcripts and morphology image for all crops in a single process.
// Replaces the separate CROP_TRANSCRIPTS + CROP_IMAGE modules.
//
// Image crops use SpatialData Image2DModel + bounding_box_query, so the
// pixel-to-micron transform is handled formally (Scale transform, pixel_size_um
// per pixel in y and x).  Transcript crops use a pandas bounding-box filter
// to preserve every column.
process CROP_SPATIALDATA {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.4'

    input:
    tuple val(meta), path(transcripts), path(morphology_tif), path(crops_csv)
    val pixel_size_um

    output:
    tuple val(meta), path("*_transcripts.parquet"), emit: cropped_transcripts
    tuple val(meta), path("*_morphology.tif"),      emit: cropped_images

    script:
    """
    crop_spatialdata.py \\
        ${transcripts} \\
        ${morphology_tif} \\
        ${crops_csv} \\
        ${pixel_size_um}
    """
    stub:
    """
    touch stub_crop1_transcripts.parquet stub_crop1_morphology.tif
    """
}
