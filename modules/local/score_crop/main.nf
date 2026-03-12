process SCORE_CROP {
    tag "${meta.id}:${crop_id}:${method}:${param_hash}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.0'

    input:
    tuple val(meta), val(crop_id), val(method), val(param_hash),
          path(cells), path(transcripts_assigned)
    val crop_area_mm2

    output:
    path "*_score.json", emit: score

    script:
    def out_name = "${meta.id}_${crop_id}_${method}_${param_hash}_score.json"
    """
    score_crop.py \\
        ${transcripts_assigned} \\
        ${cells} \\
        ${method} \\
        ${param_hash} \\
        ${crop_id} \\
        ${crop_area_mm2}
    mv score.json ${out_name}
    """
    stub:
    """
    touch ${meta.id}_${crop_id}_${method}_${param_hash}_score.json
    """
}
