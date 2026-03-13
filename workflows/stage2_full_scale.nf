// ─────────────────────────────────────────────────────────────────────────────
// Stage 2: Full-scale segmentation at optimal parameters + scoring + report
// ─────────────────────────────────────────────────────────────────────────────

include { PROSEG_FULL     } from '../modules/local/proseg_full/main'
include { CELLPOSE_FULL   } from '../modules/local/cellpose_full/main'
include { SEGGER_FULL     } from '../modules/local/segger_full/main'
include { XR_FULL         } from '../modules/local/xr_full/main'
include { SCORE_FULL      } from '../modules/local/score_full/main'
include { GENERATE_REPORT } from '../modules/local/generate_report/main'

workflow STAGE2_FULL_SCALE {
    take:
    ch_samples                // [meta, transcripts, nucleus_boundaries, morphology_tif, xenium_bundle, h5ad, xoa_version]
    ch_optimal_params         // path: optimal_params.json
    ch_scores_summary         // path: scores_summary.csv
    markers_yaml              // path: markers.yaml
    segger_model              // path or null
    mecr_weight
    recovery_weight
    yield_weight
    pixel_size_um
    nucleus_segmentation_only // bool

    main:

    // ── Parse optimal_params.json into per-method value channels ─────────────
    // We use a single-element channel from the JSON file, slurp it in Groovy.
    ch_optimal_map = ch_optimal_params.map { f ->
        new groovy.json.JsonSlurper().parse(f)
    }

    // Each method channel emits one element (the params map) if the method exists.
    ch_proseg_opt = ch_optimal_map
        .filter { m -> m.containsKey("proseg") }
        .map { m ->
            def p = m.proseg.params
            tuple(
                p.cell_compactness                as Float,
                p.max_transcript_nucleus_distance as Float,
                p.voxel_size                      as Float,
                p.diffusion_probability           as Float
            )
        }

    ch_cellpose_opt = ch_optimal_map
        .filter { m -> m.containsKey("cellpose") }
        .map { m ->
            def p = m.cellpose.params
            tuple(p.diameter as Float, p.flow_threshold as Float, p.sharpen_tiff ?: 'false')
        }

    ch_segger_opt = ch_optimal_map
        .filter { m -> m.containsKey("segger") }
        .map { m ->
            def p = m.segger.params
            // Order matches SEGGER_FULL input: min_transcripts_per_cell, tile_size
            // (module reorders to tile, min_tx internally)
            tuple(p.min_transcripts_per_cell as Integer, p.tile_size as Integer)
        }

    ch_xr_opt = ch_optimal_map
        .filter { m -> m.containsKey("xr") }
        .map { m ->
            def p = m.xr.params
            tuple(p.expansion_distance as Float, p.dapi_filter as Float, p.boundary_stain ?: 'false')
        }

    // ── Split sample channel by component ────────────────────────────────────
    ch_meta_tx     = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, tx) }
    ch_meta_nuc    = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, nuc) }
    ch_meta_img    = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, img) }
    ch_meta_bundle = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, bun) }
    ch_meta_h5ad   = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, h5) }
    ch_meta_ver    = ch_samples.map { meta, tx, nuc, img, bun, h5, ver -> tuple(meta, ver) }

    // ── ProSeg full-scale ─────────────────────────────────────────────────────
    PROSEG_FULL(
        ch_meta_tx.combine(ch_proseg_opt)
            .map { meta, tx, compact, max_nuc_dist, vox, diff ->
                tuple(meta, tx, compact, max_nuc_dist, vox, diff)
            }
    )

    // ── Cellpose full-scale ───────────────────────────────────────────────────
    CELLPOSE_FULL(
        ch_meta_tx.join(ch_meta_img).combine(ch_cellpose_opt)
            .map { meta, tx, img, diam, flow, sharpen ->
                tuple(meta, tx, img, diam, flow, sharpen)
            },
        pixel_size_um,
        nucleus_segmentation_only
    )

    // ── SEGGER full-scale (optional) ──────────────────────────────────────────
    ch_segger_results = Channel.empty()
    if (segger_model) {
        SEGGER_FULL(
            ch_meta_bundle.combine(ch_segger_opt)
                .map { meta, bundle, min_tx, tile -> tuple(meta, bundle, tile, min_tx) },
            file(segger_model)
        )
        ch_segger_results = SEGGER_FULL.out.results
    }

    // ── XR full-scale (XOA3 samples only) ────────────────────────────────────
    ch_xr3_bundles = ch_meta_bundle
        .join(ch_meta_ver, by: 0)
        .filter { meta, bundle, ver -> ver.toString() == '3' }
        .map   { meta, bundle, _ver -> tuple(meta, bundle) }

    XR_FULL(
        ch_xr3_bundles.combine(ch_xr_opt)
            .map { meta, bundle, exp, dapi, bdry ->
                tuple(meta, bundle, exp, dapi, bdry)
            },
        nucleus_segmentation_only
    )

    // ── Baseline scoring (XOA outputs from nf-xenium-processing) ─────────────
    // nucleus_boundaries.parquet is used as the "cells" file for cell count.
    // transcripts.parquet from nf-xenium-processing carries the cell_id column.
    // Method name = "xoa3" or "xoa4" based on the detected version.
    ch_baseline = ch_meta_tx
        .join(ch_meta_nuc,  by: 0)
        .join(ch_meta_h5ad, by: 0)
        .join(ch_meta_ver,  by: 0)
        .map { meta, tx, nuc, h5ad, ver ->
            tuple(meta, "xoa${ver}", nuc, tx, h5ad)
        }

    // ── Assemble all method results for scoring ───────────────────────────────
    // Format: [meta, method, cells_file, transcripts_assigned, h5ad]
    ch_proseg_s2 = PROSEG_FULL.out.results
        .join(ch_meta_h5ad, by: 0)
        .map { meta, method, cells, tx, h5ad -> tuple(meta, method, cells, tx, h5ad) }

    ch_cellpose_s2 = CELLPOSE_FULL.out.results
        .join(ch_meta_h5ad, by: 0)
        .map { meta, method, cells, tx, h5ad -> tuple(meta, method, cells, tx, h5ad) }

    ch_segger_s2 = ch_segger_results
        .join(ch_meta_h5ad, by: 0)
        .map { meta, method, cells, tx, h5ad -> tuple(meta, method, cells, tx, h5ad) }

    ch_xr_s2 = XR_FULL.out.results
        .join(ch_meta_h5ad, by: 0)
        .map { meta, method, cells, tx, h5ad -> tuple(meta, method, cells, tx, h5ad) }

    ch_all_s2 = ch_proseg_s2
        .mix(ch_cellpose_s2)
        .mix(ch_segger_s2)
        .mix(ch_xr_s2)
        .mix(ch_baseline)

    SCORE_FULL(
        ch_all_s2,
        file(markers_yaml),
        mecr_weight,
        recovery_weight,
        yield_weight,
        0      // baseline_cell_count: 0 disables yield normalization
    )

    // ── Generate HTML report ──────────────────────────────────────────────────
    GENERATE_REPORT(
        SCORE_FULL.out.score_csv.collect(),
        ch_scores_summary,
        ch_optimal_params
    )

    emit:
    report     = GENERATE_REPORT.out.report
    score_csvs = SCORE_FULL.out.score_csv
}
