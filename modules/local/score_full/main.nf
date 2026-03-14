process SCORE_FULL {
    tag "${meta.id}:${method}"
    label 'process_medium'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.4'

    publishDir "${params.outdir}/${meta.id}/scores", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(method), path(cells), path(transcripts_assigned), path(h5ad)
    path markers_yaml
    val mecr_weight
    val recovery_weight
    val yield_weight
    val baseline_cell_count  // pass 0 to disable yield normalization

    output:
    path "${meta.id}_${method}_score.csv", emit: score_csv

    script:
    def baseline_arg = (baseline_cell_count?.toString()?.toInteger() ?: 0) > 0
        ? "--baseline-cell-count ${baseline_cell_count}" : ""
    """
    score_full.py \\
        --h5ad        ${h5ad} \\
        --transcripts ${transcripts_assigned} \\
        --cells       ${cells} \\
        --markers     ${markers_yaml} \\
        --method      ${method} \\
        --sample      ${meta.id} \\
        --mecr-weight     ${mecr_weight} \\
        --recovery-weight ${recovery_weight} \\
        --yield-weight    ${yield_weight} \\
        ${baseline_arg}
    """
    stub:
    """
    touch ${meta.id}_${method}_score.csv
    """
}
