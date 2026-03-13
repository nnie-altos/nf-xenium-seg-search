// XeniumRanger resegment — Stage 2 full-slide at optimal parameters.
// Only runs when XOA version == 3.
process XR_FULL {
    tag "${meta.id}"
    label 'process_high'

    container 'nf-core/xeniumranger:4.0'

    publishDir "${params.outdir}/${meta.id}/xr", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(xenium_bundle),
          val(expansion_distance), val(dapi_filter), val(boundary_stain)
    val nucleus_segmentation_only

    output:
    tuple val(meta), val("xr"),
          path("final_output/cell_boundaries.parquet"),
          path("final_output/transcripts.parquet"),
          emit: results

    script:
    def boundary_stain_flag = boundary_stain.toString() == 'false' ? '--boundary-stain=disable' : ''
    """
    # ── Step 1: resegment (always) ───────────────────────────────────────────
    xeniumranger resegment \\
        --id            xr_output \\
        --xenium-bundle ${xenium_bundle} \\
        --expansion-distance ${expansion_distance} \\
        --dapi-filter        ${dapi_filter} \\
        ${boundary_stain_flag} \\
        --localcores    ${task.cpus} \\
        --localmem      ${task.memory.toGiga()}

    # ── Step 2: nucleus-only mode — import-segmentation ──────────────────────
    if [ "${nucleus_segmentation_only}" = "true" ]; then
        xeniumranger import-segmentation \\
            --id            xr_import_output \\
            --xenium-bundle ${xenium_bundle} \\
            --nuclei        xr_output/outs/cells.zarr.zip \\
            --expansion-distance ${expansion_distance} \\
            --units         pixels \\
            --localcores    ${task.cpus} \\
            --localmem      ${task.memory.toGiga()}
        mkdir -p final_output
        cp xr_import_output/outs/cell_boundaries.parquet final_output/
        cp xr_import_output/outs/transcripts.parquet     final_output/
    else
        mkdir -p final_output
        cp xr_output/outs/cell_boundaries.parquet final_output/
        cp xr_output/outs/transcripts.parquet     final_output/
    fi
    """
    stub:
    """
    mkdir -p final_output
    touch final_output/cell_boundaries.parquet final_output/transcripts.parquet
    """
}
