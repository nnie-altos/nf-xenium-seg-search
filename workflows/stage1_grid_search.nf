// ─────────────────────────────────────────────────────────────────────────────
// Stage 1: Parameter grid search on density-stratified crops
// ─────────────────────────────────────────────────────────────────────────────

include { GENERATE_PARAM_COMBOS  } from '../modules/local/generate_param_combos/main'
include { SELECT_CROPS           } from '../modules/local/select_crops/main'
include { CROP_SPATIALDATA       } from '../modules/local/crop_spatialdata/main'
include { PROSEG_CROP            } from '../modules/local/proseg_crop/main'
include { CELLPOSE_CROP          } from '../modules/local/cellpose_crop/main'
include { CELLPOSE_MASK_CROP     } from '../modules/local/cellpose_mask_crop/main'
include { BAYSOR_CROP            } from '../modules/local/baysor_crop/main'
include { SEGGER_PREDICT_CROP    } from '../modules/local/segger_predict_crop/main'
include { XR_GRIDSEARCH          } from '../modules/local/xr_gridsearch/main'
include { SCORE_CROP             } from '../modules/local/score_crop/main'
include { SELECT_OPTIMAL_PARAMS  } from '../modules/local/select_optimal_params/main'
include { RENDER_CROP_OVERLAY    } from '../modules/local/render_crop_overlay/main'
include { GENERATE_GRID_REPORT   } from '../modules/local/generate_grid_report/main'

workflow STAGE1_GRID_SEARCH {
    take:
    ch_samples                // [meta, transcripts, nucleus_boundaries, morphology_tif, xenium_bundle, xoa_version]
    ch_param_grids            // path: param_grids.yaml
    segger_model              // path or null
    n_crops                   // int
    crop_size_um              // float (µm)
    pixel_size_um             // float
    search_strategy           // string: 'grid' or 'coordinate_descent'
    nucleus_segmentation_only // bool

    main:

    // ── Generate param combinations ──────────────────────────────────────────
    GENERATE_PARAM_COMBOS(ch_param_grids, search_strategy, nucleus_segmentation_only)

    // ── Parse per-method param TSVs into channels of param-tuples ────────────
    def tsvs = GENERATE_PARAM_COMBOS.out.param_tsvs.flatten()

    ch_proseg_params = tsvs
        .filter { it.name.startsWith("proseg") }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.param_hash,
            row.cell_compactness                  as Float,
            row.max_transcript_nucleus_distance   as Float,
            row.voxel_size                        as Float,
            row.diffusion_probability             as Float
        )}

    ch_cellpose_params = tsvs
        .filter { it.name.startsWith("cellpose") }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.param_hash,
            row.diameter        as Float,
            row.flow_threshold  as Float,
            row.sharpen_tiff ?: 'false'
        )}

    ch_segger_params = tsvs
        .filter { it.name.startsWith("segger") }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.param_hash,
            row.min_transcripts_per_cell as Integer,
            row.dist_tx                  as Float,
            row.tile_size                as Integer
        )}

    ch_xr_params = tsvs
        .filter { it.name.startsWith("xr") }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.param_hash,
            row.expansion_distance as Float,
            row.dapi_filter        as Float,
            row.boundary_stain ?: 'false'
        )}

    ch_baysor_params = tsvs
        .filter { it.name.startsWith("cellpose_baysor") }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.param_hash,
            row.prior_segmentation_confidence as Float,
            row.min_molecules_per_cell        as Integer
        )}

    // ── Split sample channel by component ────────────────────────────────────
    ch_meta_tx      = ch_samples.map { meta, tx, nuc, img, bun, ver -> tuple(meta, tx) }
    ch_meta_nucleus = ch_samples.map { meta, tx, nuc, img, bun, ver -> tuple(meta, nuc) }
    ch_meta_image   = ch_samples.map { meta, tx, nuc, img, bun, ver -> tuple(meta, img) }
    ch_meta_bundle  = ch_samples.map { meta, tx, nuc, img, bun, ver -> tuple(meta, bun) }
    ch_meta_ver     = ch_samples.map { meta, tx, nuc, img, bun, ver -> tuple(meta, ver) }

    // ── Select density-stratified crops ──────────────────────────────────────
    SELECT_CROPS(ch_meta_nucleus, n_crops, crop_size_um, pixel_size_um)
    ch_crops = SELECT_CROPS.out.crops  // [meta, crops.csv]

    // ── Crop transcripts + images (SpatialData bounding_box_query) ───────────
    CROP_SPATIALDATA(
        ch_meta_tx
            .join(ch_meta_image, by: 0)
            .join(ch_crops,      by: 0)
            .map { meta, tx, img, crops -> tuple(meta, tx, img, crops) },
        pixel_size_um
    )
    // Flatten lists of per-crop files to individual [meta, crop_id, file] tuples
    ch_cropped_tx = CROP_SPATIALDATA.out.cropped_transcripts
        .transpose()
        .map { meta, tx ->
            def crop_id = tx.name.replaceFirst(/_transcripts\.parquet$/, '')
            tuple(meta, crop_id, tx)
        }
    ch_cropped_img = CROP_SPATIALDATA.out.cropped_images
        .transpose()
        .map { meta, img ->
            def crop_id = img.name.replaceFirst(/_morphology\.tif$/, '')
            tuple(meta, crop_id, img)
        }

    // ── Compute crop area in mm² ──────────────────────────────────────────────
    def crop_size_mm  = crop_size_um.toFloat() / 1000.0
    def crop_area_mm2 = crop_size_mm * crop_size_mm

    // ── ProSeg crop grid search ───────────────────────────────────────────────
    PROSEG_CROP(
        ch_cropped_tx.combine(ch_proseg_params)
            .map { meta, crop_id, tx, ph, compact, max_nuc_dist, vox, diff ->
                tuple(meta, crop_id, tx, ph, compact, max_nuc_dist, vox, diff)
            }
    )

    // ── Cellpose crop grid search ─────────────────────────────────────────────
    // Join tx + img + crops_csv by (meta, crop_id) then cross with params.
    // crops.csv is joined by meta so each sample uses its own crops.csv.
    ch_tx_img_csv = ch_cropped_tx
        .join(ch_cropped_img, by: [0, 1])           // [meta, crop_id, tx, img]
        .join(ch_crops.map { meta, csv -> [meta, csv] }, by: 0)  // add crops.csv keyed by meta
        // → [meta, crop_id, tx, img, crops.csv]

    CELLPOSE_CROP(
        ch_tx_img_csv.combine(ch_cellpose_params)
            .map { meta, crop_id, tx, img, crops_csv, ph, diam, flow, sharpen ->
                tuple(meta, crop_id, tx, img, crops_csv, ph, diam, flow, sharpen)
            },
        pixel_size_um,
        nucleus_segmentation_only
    )

    // ── Cellpose+Baysor crop grid search ─────────────────────────────────────
    // Cellpose runs once per crop at default params to generate the prior mask.
    // Only Baysor params (prior_segmentation_confidence, min_molecules_per_cell) are varied.
    // Default diameter: nucleus mode = 15 px, cell mode = 30 px.
    def baysor_default_diameter = nucleus_segmentation_only.toString() == 'true' ? 15.0f : 30.0f

    CELLPOSE_MASK_CROP(
        ch_cropped_img,
        nucleus_segmentation_only,
        baysor_default_diameter,
        0.4   // default flow_threshold
    )

    ch_baysor_crop_inputs = CELLPOSE_MASK_CROP.out.mask
        .join(ch_cropped_tx, by: [0, 1])
        .map { meta, crop_id, mask, tx -> tuple(meta, crop_id, tx, mask) }
        .join(ch_crops.map { meta, csv -> [meta, csv] }, by: 0)
        .map { meta, crop_id, tx, mask, crops_csv -> tuple(meta, crop_id, tx, mask, crops_csv) }
        .combine(ch_baysor_params)
        .map { meta, crop_id, tx, mask, crops_csv, ph, psc, mmpc ->
            tuple(meta, crop_id, tx, mask, crops_csv, ph, psc, mmpc)
        }

    BAYSOR_CROP(ch_baysor_crop_inputs, pixel_size_um, baysor_default_diameter)

    // ── SEGGER grid search (optional — only if model provided) ───────────────
    // SEGGER runs on the full Xenium bundle (not crops): it needs
    // nucleus_boundaries.parquet for spatial graph construction.
    // tile_size → dataset creation; min_transcripts_per_cell → post-processing.
    ch_segger_results = Channel.empty()
    if (segger_model) {
        SEGGER_PREDICT_CROP(
            ch_meta_bundle.combine(ch_segger_params)
                .map { meta, bundle, ph, min_tx, dist_tx, tile ->
                    tuple(meta, bundle, ph, tile, dist_tx, min_tx)
                },
            file(segger_model)
        )
        ch_segger_results = SEGGER_PREDICT_CROP.out.results
    }

    // ── XeniumRanger grid search (XOA3 only — full-slide) ────────────────────
    ch_xr3_bundles = ch_meta_bundle
        .join(ch_meta_ver, by: 0)
        .filter { meta, bundle, ver -> ver.toString() == '3' }
        .map   { meta, bundle, _ver -> tuple(meta, bundle) }

    XR_GRIDSEARCH(
        ch_xr3_bundles.combine(ch_xr_params)
            .map { meta, bundle, ph, exp, dapi, bdry ->
                tuple(meta, bundle, ph, exp, dapi, bdry)
            }
    )

    // ── Mix all segmentation results and score ────────────────────────────────
    // All results share the same tuple shape:
    // [meta, crop_id, method, param_hash, cells_file, transcripts_assigned]
    ch_all_seg = PROSEG_CROP.out.results
        .mix(CELLPOSE_CROP.out.results)
        .mix(BAYSOR_CROP.out.results)
        .mix(ch_segger_results)
        .mix(XR_GRIDSEARCH.out.results)

    SCORE_CROP(ch_all_seg, crop_area_mm2)

    // ── Collect all scores and select optimal params ──────────────────────────
    SELECT_OPTIMAL_PARAMS(
        SCORE_CROP.out.score.collect(),
        GENERATE_PARAM_COMBOS.out.manifest
    )

    // ── Render per-combo overlays (DAPI + cell boundaries coloured by tx count) ─
    // Join segmentation results with their crop images by (meta, crop_id)
    ch_overlay_inputs = ch_all_seg
        .join(ch_cropped_img, by: [0, 1])
        // → [meta, crop_id, method, param_hash, cells, tx, morphology_tif]

    RENDER_CROP_OVERLAY(ch_overlay_inputs)

    // ── Assemble grid search HTML report ─────────────────────────────────────
    GENERATE_GRID_REPORT(
        SCORE_CROP.out.score.collect(),
        RENDER_CROP_OVERLAY.out.png.collect()
    )

    emit:
    optimal_params    = SELECT_OPTIMAL_PARAMS.out.optimal_params
    scores_summary    = SELECT_OPTIMAL_PARAMS.out.scores_summary
    crops             = ch_crops
    grid_report       = GENERATE_GRID_REPORT.out.report
}
