// Pairwise IoU-based Average Precision scoring across all Stage 2 methods.
// Input channel must already be grouped by sample (groupTuple).
// Format: [meta, [method1, method2, ...], [cells1.parquet, cells2.parquet, ...]]
process SCORE_AP {
    tag "${meta.id}"
    label 'process_medium'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.2'

    publishDir "${params.outdir}/${meta.id}/scores", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(methods), path(cells_files)

    output:
    path "${meta.id}_ap_matrix.csv", emit: ap_matrix

    script:
    def method_args = [methods, cells_files instanceof List ? cells_files : [cells_files]]
        .transpose()
        .collect { m, f -> "${m}=${f}" }
        .join(" ")
    """
    score_ap.py ${meta.id} ${method_args}
    """
    stub:
    """
    printf 'sample_id,method_a,method_b,n_cells_a,n_cells_b,ap,ap_at_0.5,ap_at_0.8\\n' > ${meta.id}_ap_matrix.csv
    """
}
