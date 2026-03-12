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

    output:
    tuple val(meta), val("xr"),
          path("xr_output/outs/cell_boundaries.parquet"),
          path("xr_output/outs/transcripts.parquet"),
          emit: results

    script:
    // Default behaviour includes boundary stain; explicitly disable when not wanted
    def boundary_stain_flag = boundary_stain.toString() == 'false' ? '--boundary-stain=disable' : ''
    """
    xeniumranger resegment \\
        --id            xr_output \\
        --xenium-bundle ${xenium_bundle} \\
        --expansion-distance ${expansion_distance} \\
        --dapi-filter        ${dapi_filter} \\
        ${boundary_stain_flag} \\
        --localcores    ${task.cpus} \\
        --localmem      ${task.memory.toGiga()}
    """
    stub:
    """
    mkdir -p xr_output/outs && touch xr_output/outs/cell_boundaries.parquet xr_output/outs/transcripts.parquet
    """
}
