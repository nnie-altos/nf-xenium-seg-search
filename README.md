# nf-xenium-seg-search

A Nextflow pipeline for automated cell segmentation parameter optimisation and cross-method comparison on 10x Genomics Xenium spatial transcriptomics data.

## Overview

Takes the outputs of [nf-xenium-processing](https://github.com/altos-labs/nf-xenium-processing) as input and performs:

**Stage 1 — Parameter grid search** on density-stratified spatial crops:
- Selects `n_crops` × `crop_size_um` µm crops proportionally from low / mid / high density regions
- Runs ProSeg, Cellpose, and SEGGER (optional) with all parameter combinations (full grid or coordinate descent)
- XeniumRanger resegment (XOA3 only) runs as a fixed baseline — not grid-searched
- Scores each combination by: transcript assignment rate (60%) + normalised cell yield (40%); reports mean transcripts/cell as a diagnostic
- Selects the best parameter set per method and generates an HTML grid search report with DAPI + boundary overlays

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

CSV with three columns (all other paths are derived automatically):

| Column | Description |
|--------|-------------|
| `id` | Unique sample identifier |
| `xenium_bundle` | Path to the Xenium output bundle directory |
| `nf_outdir` | Output directory of the upstream nf-xenium-processing run |

Derived paths (resolved at runtime):
- `transcripts.parquet`, `nucleus_boundaries.parquet`, `experiment.xenium` — from `xenium_bundle`
- `morphology_focus_0000.ome.tif` (XOA3) or `ch0000_dapi.ome.tif` (XOA4) — detected automatically
- `spatial_with_annotations.h5ad` — from `${nf_outdir}/${id}/${id}/`

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
| `--search_strategy` | `grid` | Parameter search strategy: `grid` (full factorial) or `coordinate_descent` |
| `--nucleus_segmentation_only` | `false` | Segment nuclei only (affects Cellpose model/diameter and XR import-segmentation) |
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
| Python scripts | `ghcr.io/nnie-altos/xenium-seg-search:0.1.0` (see `Dockerfile`) |
| ProSeg | `ghcr.io/dcjones/proseg:v3.1.0` |
| Cellpose | `ghcr.io/mouseland/cellpose:3.0.11` |
| SEGGER | `danielunyi42/segger_cuda118_py311:latest` |
| XeniumRanger | `nf-core/xeniumranger:4.0` |

Build the custom container:
```bash
docker build -t ghcr.io/nnie-altos/xenium-seg-search:0.1.0 .
```

## Output structure

```
results/
├── pipeline_info/                  # Nextflow execution reports
├── stage1/
│   ├── optimal_params.json         # Best params per method
│   ├── scores_summary.csv          # All Stage 1 scores
│   └── grid_search_report.html     # Visual report: top/bottom combos per method with DAPI overlays
├── {sample_id}/
│   ├── crops/
│   │   └── crops.csv
│   ├── proseg/
│   ├── cellpose/
│   ├── segger/                     # if --segger_model provided
│   ├── xr/                         # if XOA3 sample
│   └── scores/
│       └── *_score.csv
└── segmentation_report.html
```
