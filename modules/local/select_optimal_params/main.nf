process SELECT_OPTIMAL_PARAMS {
    tag "select_optimal_params"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.4'

    publishDir "${params.outdir}/stage1", mode: params.publish_dir_mode

    input:
    path score_jsons  // collected list of all score.json files
    path param_manifest

    output:
    path "optimal_params.json", emit: optimal_params
    path "scores_summary.csv",  emit: scores_summary

    script:
    """
    aggregate_scores.py ${score_jsons} --param-manifest ${param_manifest}
    """
    stub:
    """
    echo '{}' > optimal_params.json
touch scores_summary.csv
    """
}
