// Assembles the Stage 1 grid search HTML report.
// Selects top-3 and bottom-3 parameter combinations per method per crop
// based on composite score, and embeds pre-rendered overlay PNGs.
process GENERATE_GRID_REPORT {
    tag "grid_report"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.3'

    publishDir "${params.outdir}/stage1", mode: params.publish_dir_mode

    input:
    path score_jsons   // collected list of all *_score.json files
    path overlay_pngs  // collected list of all *_overlay.png files

    output:
    path "grid_search_report.html", emit: report

    script:
    """
    generate_grid_report.py ${score_jsons} --overlays ${overlay_pngs}
    """
    stub:
    """
    touch grid_search_report.html
    """
}
