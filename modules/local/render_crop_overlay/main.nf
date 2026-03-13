// Renders a DAPI crop with cell boundary outlines coloured by transcript count.
// One process invocation per (method × param_hash × crop_id).
// Output PNG filename encodes all metadata for downstream report assembly.
process RENDER_CROP_OVERLAY {
    tag "${meta.id}:${crop_id}:${method}:${param_hash}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.1'

    input:
    tuple val(meta), val(crop_id), val(method), val(param_hash),
          path(cells), path(transcripts_assigned), path(morphology_tif)

    output:
    path "*_overlay.png", emit: png

    script:
    """
    render_crop_overlay.py \\
        ${morphology_tif} \\
        ${cells} \\
        ${transcripts_assigned} \\
        ${method} \\
        ${param_hash} \\
        ${crop_id} \\
        ${meta.id}
    """
    stub:
    """
    touch ${meta.id}_${crop_id}_${method}_${param_hash}_overlay.png
    """
}
