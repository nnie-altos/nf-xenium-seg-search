// XeniumRanger resegment — Stage 1 grid search.
// XR does not support spatial cropping, so each combination runs on the full
// Xenium bundle. Scoring is handled by SCORE_CROP (same process) but the
// "crop_area_mm2" is set to the full slide area and crop_id = "fullslide".
// Only runs when XOA version == 3.
process XR_GRIDSEARCH {
    tag "${meta.id}:${param_hash}"
    label 'process_high'

    container 'nf-core/xeniumranger:4.0'

    input:
    tuple val(meta), path(xenium_bundle), val(param_hash),
          val(expansion_distance), val(dapi_filter), val(boundary_stain)

    output:
    tuple val(meta), val("fullslide"), val("xr"), val(param_hash),
          path("xr_output/outs/cell_boundaries.parquet"),
          path("xr_output/outs/transcripts.parquet"),
          emit: results

    script:
    // Default behaviour includes boundary stain; explicitly disable when not wanted
    def boundary_stain_flag = boundary_stain.toString() == 'false' ? '--boundary-stain=disable' : ''
    """
    xeniumranger resegment \\
        --id           xr_output \\
        --xenium-bundle ${xenium_bundle} \\
        --expansion-distance ${expansion_distance} \\
        --dapi-filter        ${dapi_filter} \\
        ${boundary_stain_flag} \\
        --localcores   ${task.cpus} \\
        --localmem     ${task.memory.toGiga()}
    """
}
