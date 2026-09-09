# S1 · Data and preprocessing

Aligns nine snATAC-seq samples with Cell Ranger, filters barcodes, removes
doublets, and produces the annotated `AnnDataSet` used for pseudobulk peak
calling in S2.

## Scripts (run in order)

| # | File | Purpose |
|---|------|---------|
| 01 | `01_cellranger_atac.sh` | 9-sample `cellranger-atac count` loop |
| 02 | `02_qc_filter.ipynb`    | `snap.pp.filter_cells` — TSSe ≥ 6, `passed_filters + 1 ∈ [5000, 50000]`, `is__cell_barcode == 1` |
| 03 | `03_doublet_removal.ipynb` | `snap.pp.add_tile_matrix` (500 bp) → `snap.pp.select_features(n_features=250000)` → `snap.pp.scrublet` → `snap.pp.filter_doublets(probability_threshold=0.5)` |
| 04 | `04_clustering_annotation.ipynb` | `snap.pp.select_features(n_features=200000)` → `snap.tl.spectral(30)` → `snap.pp.harmony(batch="sample", max_iter_harmony=20)` → `snap.pp.knn(k=10)` → `snap.tl.leiden(random_state=0)` → manual cell-type labels |

## Cell Ranger provenance

Binary: `cellranger-atac 2.1.0` (ATAC-only pipeline)
Chemistry flag: `--chemistry=ARC-v1` (multiome ATAC library format)
Reference: `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` (ARC reference used
with the ATAC binary; `mkref_version=cellranger-arc-2.0.0` per each sample's
`fragments.tsv.gz` header).

Verified from all 9 samples' `_cmdline`, `_versions`, `outs/summary.csv`, and
`outs/fragments.tsv.gz` headers — every sample used identical binary /
version / reference / chemistry / hashes.

## `n_features` — 250000 vs 200000 (intentional divergence)

`03_doublet_removal.ipynb` runs `select_features(n_features=250000)` per
sample (larger feature set for Scrublet's per-sample simulation), while
`04_clustering_annotation.ipynb` uses `select_features(n_features=200000)`
on the merged `AnnDataSet` before spectral embedding. Both values are the
on-disk numbers and are preserved as-is.

## Inputs

- Raw fastqs at `../data/sample{1..9}/*.fastq.gz`
- `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` reference bundle

## Outputs (not tracked; see `data/README.md`)

- `cellranger_output_heart/sample{1..9}/outs/` (per Cell Ranger structure)
- Post-QC AnnData per sample after step 02
- Post-doublet AnnData per sample after step 03
- Merged, annotated `AnnDataSet` after step 04 (44,262 nuclei, 10 cell types)

## Environment

`envs/snapatac2.yml` (SnapATAC2 2.8.0, harmonypy 0.0.10, scanpy 1.10.3;
Python 3.9.21).
