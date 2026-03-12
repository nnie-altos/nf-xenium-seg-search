# nf-xenium-seg-search

A Nextflow pipeline for automated cell segmentation parameter optimisation and cross-method comparison on 10x Genomics Xenium spatial transcriptomics data.

## Overview

Takes the outputs of [nf-xenium-processing](https://github.com/altos-labs/nf-xenium-processing) as input and performs:

**Stage 1 — Parameter grid search** on density-stratified spatial crops:
- Selects 10 × 500 µm crops proportionally from low / mid / high density regions
- Runs ProSeg, Cellpose, SEGGER (optional), and XeniumRanger resegment (XOA3 only) with all parameter combinations
- Scores each combination by: transcript assignment rate (60%) + normalised cell yield (40%)
- Selects the best parameter set per method

**Stage 2 — Full-slide segmentation** at optimal parameters:
- Runs all methods at their best parameters on the full slide
- Scores by: MECR (50%) + transcript recovery (30%) + cell yield (20%)
- Generates an HTML report with comparative plots and ROI overlays

## Requirements

- Nextflow ≥ 23.04
- Docker (or Singularity)
- GPU recommended for Cellpose and SEGGER

## Quick start

```bash
nextflow run main.nf \
    --input  assets/test/samplesheet.csv \
    --outdir results/ \
    --markers assets/markers/lung_markers.yaml \
    -profile docker,test
```

## Samplesheet format

CSV with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `transcripts` | Path to `transcripts.parquet` (nf-xenium-processing output) |
| `nucleus_boundaries` | Path to `nucleus_boundaries.parquet` (from Xenium bundle) |
| `morphology_tif` | Path to `morphology_focus.ome.tif` (from Xenium bundle) |
| `xenium_bundle` | Path to full Xenium bundle directory |
| `experiment_xenium` | Path to `experiment.xenium` (XOA version detection) |
| `h5ad` | Path to `spatial_with_annotations.h5ad` (nf-xenium-processing output) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Samplesheet CSV (required) |
| `--outdir` | — | Output directory (required) |
| `--markers` | — | Path to `markers.yaml` (required) |
| `--segger_model` | `null` | Path to pre-trained SEGGER model (SEGGER skipped if absent) |
| `--crop_size_um` | `500` | Crop side length in microns |
| `--n_crops` | `10` | Number of density-stratified crops |
| `--pixel_size_um` | `0.2125` | Xenium pixel size at Level 0 (µm/pixel) |
| `--mecr_weight` | `0.5` | MECR weight in Stage 2 composite score |
| `--recovery_weight` | `0.3` | Transcript recovery weight |
| `--yield_weight` | `0.2` | Cell yield weight |
| `--skip_stage1` | `false` | Skip parameter search; use `--optimal_params` |
| `--optimal_params` | `null` | Path to `optimal_params.json` from a prior Stage 1 run |
| `--param_grids` | `conf/param_grids.yaml` | Parameter grid definitions |

## Markers file format

```yaml
markers:
  GeneA: CellTypeX   # gene → cell type it marks
  GeneB: CellTypeX
  GeneC: CellTypeY
  GeneD: CellTypeZ
```

Genes are fuzzy-matched against the h5ad `var_names` (case-insensitive, ignoring `-_` and spaces). See `assets/markers/lung_markers.yaml` for an example.

## Containers

| Process | Container |
|---------|-----------|
| Python scripts | `ghcr.io/altos-labs/xenium-seg-search:0.1.0` (see `Dockerfile`) |
| ProSeg | `ghcr.io/dcjones/proseg:v3.1.0` |
| Cellpose | `ghcr.io/mouseland/cellpose:3.0.11` |
| SEGGER | `danielunyi42/segger_cuda118_py311:latest` |
| XeniumRanger | `nfcore/xeniumranger:3.1.1` |

Build the custom container:
```bash
docker build -t ghcr.io/altos-labs/xenium-seg-search:0.1.0 .
```

## Output structure

```
results/
├── pipeline_info/           # Nextflow execution reports
├── stage1/
│   ├── optimal_params.json  # Best params per method
│   └── scores_summary.csv   # All Stage 1 scores
├── {sample_id}/
│   ├── crops/
│   │   └── crops.csv
│   ├── proseg/
│   ├── cellpose/
│   ├── segger/              # if --segger_model provided
│   ├── xr/                  # if XOA3 sample
│   └── scores/
│       └── *_score.csv
└── segmentation_report.html
```
