# S1 · Data and preprocessing

Aligns nine snATAC-seq samples with Cell Ranger, filters barcodes, removes
doublets, clusters and annotates cell types, and exports per-(sample × cell
type) fragment BEDs consumed by S2.

## Scripts (run in order)

| # | File | Purpose |
|---|------|---------|
| 01 | `01_cellranger_atac.sh` | 9-sample `cellranger-atac count` loop |
| 02 | `02_qc_filter.py` | TSSe ≥ 6, `passed_filters + 1 ∈ [5000, 50000]`, `is__cell_barcode == 1`; writes per-sample filtered CSV and cached h5ad. Optional QC figures via `--qc-figures` |
| 03 | `03_doublet_removal.py` | `snap.pp.add_tile_matrix` (500 bp) → `snap.pp.select_features(n_features=250000)` → `snap.pp.scrublet` → `snap.pp.filter_doublets` (probability > 0.5) |
| 04 | `04_clustering_annotation.py` | Merge per-sample h5ads → `select_features(n_features=200000)` → `snap.tl.spectral(30)` → `snap.pp.harmony(batch="sample", max_iter_harmony=20)` → `snap.pp.knn(k=10)` → `snap.tl.leiden(random_state=0)` → gene-matrix dotplot → cluster→celltype mapping → per-`sample_cluster` fragment BEDs |

### Config files

- `celltype_mapping.json` — Leiden cluster id (string) → cell-type name. Assigned manually from the marker-gene dotplot; unmapped clusters become `Unknown`. Edit to match your own clustering.
- `marker_genes.json` — marker gene groups used to render the dotplot in step 04.

### Example invocation (default paths from `$PROJECT_ROOT`/`$DATA_ROOT`)

```
export PROJECT_ROOT=/path/to/your/project
export DATA_ROOT=$PROJECT_ROOT/output_so

bash  01_cellranger_atac.sh
python 02_qc_filter.py --qc-figures
python 03_doublet_removal.py
python 04_clustering_annotation.py \
    --celltype-mapping celltype_mapping.json \
    --marker-genes    marker_genes.json
```

## Cell Ranger provenance

- Binary: `cellranger-atac 2.1.0` (ATAC-only pipeline)
- Chemistry flag: `--chemistry=ARC-v1` (multiome ATAC library format)
- Reference: `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` (ARC reference used with the ATAC binary; `mkref_version=cellranger-arc-2.0.0` per each sample's `fragments.tsv.gz` header)

Verified from all 9 samples' `_cmdline`, `_versions`, `outs/summary.csv`, and `outs/fragments.tsv.gz` headers — every sample used identical binary / version / reference / chemistry / hashes.

## `n_features` — 250000 vs 200000 (intentional divergence)

`03_doublet_removal.py` runs `select_features(n_features=250000)` per sample (larger feature set for Scrublet's per-sample simulation), while `04_clustering_annotation.py` uses `select_features(n_features=200000)` on the merged `AnnDataSet` before spectral embedding. Both values match the manuscript run and are preserved as-is.

## Inputs

- Raw fastqs at `$FASTQ_ROOT/sample{1..9}/*.fastq.gz`
- `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` reference bundle

## Outputs (not tracked; see `data/README.md`)

- `cell_ranger_output/sample{1..9}/outs/` (per Cell Ranger structure)
- Per-sample filtered CSV (`$DATA_ROOT/filtered_samples/sample{i}_filtered.csv`)
- Cached per-sample h5ad after QC (`$DATA_ROOT/h5ad/sample{i}.h5ad`)
- Doublet-filtered h5ad per sample (`$DATA_ROOT/cell_qc/doublets/sample{i}_doublets.h5ad`)
- Merged, annotated `AnnDataSet` after step 04 (44,262 nuclei, 10 cell types)
- Per-(sample × cell type) fragment BEDs at `$DATA_ROOT/fragments/mcluster<CELLTYPE>.<SAMPLE>_fragments.bed` (S2 input)

## Environment

`envs/snapatac2.yml` (SnapATAC2 2.8.0, harmonypy 0.0.10, scanpy 1.10.3; Python 3.9.21).
