process GENERATE_REPORT {
    tag "report"
    label 'process_medium'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.1'

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path score_csvs     // collected Stage 2 score CSVs (one per sample × method)
    path scores_summary // Stage 1 scores_summary.csv from SELECT_OPTIMAL_PARAMS
    path optimal_params // Stage 1 optimal_params.json (= param_manifest for report)

    output:
    path "seg_search_report.html", emit: report

    script:
    """
    mkdir -p scores_dir
    cp ${score_csvs} scores_dir/

    generate_report.py \\
        --scores-dir      scores_dir \\
        --stage1-summary  ${scores_summary} \\
        --stage1-manifest ${optimal_params} \\
        --outdir          .
    """
    stub:
    """
    touch seg_search_report.html
    """
}
