process CROP_TRANSCRIPTS {
    tag "${meta.id}"
    label 'process_low'

    container 'ghcr.io/nnie-altos/xenium-seg-search:0.1.1'

    input:
    tuple val(meta), path(transcripts), path(crops_csv)

    output:
    tuple val(meta), path("*_transcripts.parquet"), emit: cropped_transcripts

    script:
    """
    crop_transcripts.py ${crops_csv} ${transcripts}
    """
    stub:
    """
    touch stub_crop1_transcripts.parquet
    """
}
