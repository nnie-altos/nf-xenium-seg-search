// Export the winning segmentation method and its optimal parameters as a
// nf-xenium-processing-compatible params YAML (-params-file).
process EXPORT_PARAMS {
    tag "export_params"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.3'

    publishDir "${params.outdir}/stage1", mode: params.publish_dir_mode

    input:
    path  score_csvs              // collected Stage 2 *_score.csv files
    path  optimal_params          // Stage 1 optimal_params.json
    val   nucleus_segmentation_only

    output:
    path "nfxp_params.yaml",   emit: nfxp_params
    path "winner_summary.json", emit: winner_summary

    script:
    def nucleus_flag = nucleus_segmentation_only.toString()
    """
    mkdir -p scores_dir
    cp ${score_csvs} scores_dir/

    export_params.py \\
        --scores-dir                 scores_dir \\
        --optimal-params             ${optimal_params} \\
        --nucleus-segmentation-only  ${nucleus_flag}
    """
    stub:
    """
    echo 'segmentation: proseg' > nfxp_params.yaml
    echo '{}' > winner_summary.json
    """
}
