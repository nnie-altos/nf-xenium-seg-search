#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// nf-xenium-seg-search
// Automated parameter grid search and cross-method comparison for
// 10x Genomics Xenium cell segmentation.
//
// Usage:
//   nextflow run main.nf \
//     --input  samplesheet.csv \
//     --outdir results/ \
//     --markers assets/markers/lung_markers.yaml
//
// Samplesheet columns (CSV):
//   id, xenium_bundle, nf_outdir
//   (all other paths are derived from these two directories)
//
// See README.md and conf/param_grids.yaml for configuration.
// ─────────────────────────────────────────────────────────────────────────────

include { DETECT_XOA_VERSION   } from './modules/local/detect_version/main'
include { STAGE1_GRID_SEARCH   } from './workflows/stage1_grid_search'
include { STAGE2_FULL_SCALE    } from './workflows/stage2_full_scale'

// ── Validate required params ─────────────────────────────────────────────────
def validateParams() {
    if (!params.input)   error "ERROR: --input (samplesheet CSV) is required"
    if (!params.outdir)  error "ERROR: --outdir is required"
    if (!params.markers) error "ERROR: --markers (markers.yaml) is required"
}

// ── Parse samplesheet CSV ────────────────────────────────────────────────────
// Expected columns: id, xenium_bundle, nf_outdir
// Derived paths:
//   transcripts        = ${xenium_bundle}/transcripts.parquet
//   nucleus_boundaries = ${xenium_bundle}/nucleus_boundaries.parquet
//   morphology_tif     = ${xenium_bundle}/morphology_focus/morphology_focus_0000.ome.tif (XOA3)
//                     or ${xenium_bundle}/morphology_focus/ch0000_dapi.ome.tif       (XOA4)
//   experiment_xenium  = ${xenium_bundle}/experiment.xenium
//   h5ad               = ${nf_outdir}/${id}/${id}/spatial_with_annotations.h5ad
def parseSamplesheet(csv_path) {
    Channel.fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta    = [ id: row.id ]
            def bundle  = row.xenium_bundle.trim()
            def nf_out  = row.nf_outdir.trim()

            // Resolve morphology focus — XOA3 vs XOA4 naming differs
            def morph_v3 = file("${bundle}/morphology_focus/morphology_focus_0000.ome.tif")
            def morph_v4 = file("${bundle}/morphology_focus/ch0000_dapi.ome.tif")
            def morph    = morph_v3.exists() ? morph_v3 : morph_v4

            [
                meta,
                file("${bundle}/transcripts.parquet",        checkIfExists: true),
                file("${bundle}/nucleus_boundaries.parquet",  checkIfExists: true),
                morph,
                file(bundle,                                  checkIfExists: true),
                file("${bundle}/experiment.xenium",           checkIfExists: true),
                file("${nf_out}/${row.id}/${row.id}/spatial_with_annotations.h5ad", checkIfExists: true)
            ]
        }
}

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {
    validateParams()

    ch_samplesheet = parseSamplesheet(params.input)

    // ── Detect XOA version per sample ─────────────────────────────────────────
    ch_experiment_xenium = ch_samplesheet
        .map { meta, tx, nuc, img, bundle, exp_xen, h5ad -> [meta, exp_xen] }

    DETECT_XOA_VERSION(ch_experiment_xenium)
    // Output: [meta, xoa_version_string]  (e.g. "3" or "4")

    // ── Build combined sample channel for Stage 1 ─────────────────────────────
    // [meta, transcripts, nucleus_boundaries, morphology_tif, xenium_bundle, xoa_version]
    ch_s1_samples = ch_samplesheet
        .map { meta, tx, nuc, img, bundle, _exp, _h5ad -> [meta, tx, nuc, img, bundle] }
        .join(DETECT_XOA_VERSION.out.version, by: 0)
        .map { meta, tx, nuc, img, bundle, ver -> [meta, tx, nuc, img, bundle, ver] }

    // ── Stage 1: grid search ──────────────────────────────────────────────────
    if (!params.skip_stage1 && !params.optimal_params) {

        STAGE1_GRID_SEARCH(
            ch_s1_samples,
            file(params.param_grids),
            params.segger_model ? file(params.segger_model) : null,
            params.n_crops,
            params.crop_size_um,
            params.pixel_size_um,
            params.search_strategy,
            params.nucleus_segmentation_only
        )

        ch_optimal_params  = STAGE1_GRID_SEARCH.out.optimal_params
        ch_scores_summary  = STAGE1_GRID_SEARCH.out.scores_summary

    } else {
        // Skip Stage 1 and use supplied optimal params
        ch_optimal_params = Channel.value(
            file(params.optimal_params ?: "${params.outdir}/stage1/optimal_params.json")
        )
        ch_scores_summary = Channel.value(
            file(params.optimal_params ?
                 "${file(params.optimal_params).parent}/scores_summary.csv" :
                 "${params.outdir}/stage1/scores_summary.csv")
        )
    }

    // ── Stage 2: full-scale segmentation + scoring + report ───────────────────
    // [meta, transcripts, nucleus_boundaries, morphology_tif, xenium_bundle, h5ad, xoa_version]
    ch_s2_samples = ch_samplesheet
        .map { meta, tx, nuc, img, bundle, _exp, h5ad -> [meta, tx, nuc, img, bundle, h5ad] }
        .join(DETECT_XOA_VERSION.out.version, by: 0)
        .map { meta, tx, nuc, img, bundle, h5ad, ver -> [meta, tx, nuc, img, bundle, h5ad, ver] }

    STAGE2_FULL_SCALE(
        ch_s2_samples,
        ch_optimal_params,
        ch_scores_summary,
        params.markers,
        params.segger_model ? file(params.segger_model) : null,
        params.mecr_weight,
        params.recovery_weight,
        params.yield_weight,
        params.pixel_size_um,
        params.nucleus_segmentation_only
    )
}

// ── Workflow completion handler ───────────────────────────────────────────────
workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║        nf-xenium-seg-search completed                ║
    ╚══════════════════════════════════════════════════════╝
    Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration : ${workflow.duration}
    Output   : ${params.outdir}
    """.stripIndent()
}
