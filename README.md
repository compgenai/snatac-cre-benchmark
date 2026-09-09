# snatac-cre-benchmark

Analysis code for **"A comprehensive guide for identifying *cis*-regulatory elements in multi-sample single-cell epigenomic data"** (Cho et al.).

This repository reproduces the benchmark reported in the Benchmark section and Figure 4 of the manuscript, in which two peak-identification workflows are applied to the same single-nucleus ATAC-seq dataset (nine human embryonic heart samples from ENCODE) and compared.

- **Standard workflow** — MACS2 at the default threshold (q ≤ 0.05), merged per cell type by fixed-width iterative overlap removal (`snap.tl.merge_peaks`, 501 bp).
- **Suggested workflow** — MACS2 at a permissive threshold (p ≤ 0.01), 3D clustering over (start, summit, end), Weiszfeld geometric-median centroids with asymmetric intervals, 50% reciprocal-overlap deduplication, and gkmQC AUROC-based trimming at 0.7.

Section numbering below mirrors **Additional file 1 (Supplementary Methods)** of the manuscript.

---

## Repository layout

```
envs/                 conda environment files (four separate environments)
s1_preprocessing/     S1  Data and preprocessing
s2_peak_calling/      S2  Peak calling and construction of cell-type peak sets
s3_gkmqc/             S3  Sequence-based quality assessment and trimming
s4_peak_comparison/   S4  Comparison of peak sets
s5_ldsc/              S5  Partitioned heritability
figures/              Scripts that generate the manuscript figures
data/                 Inputs (not tracked; see data/README.md)
results/              Outputs (not tracked; see results/README.md)
```

---

## S1. Data and preprocessing

| Script | Purpose |
|---|---|
| `01_cellranger_atac.sh` | `cellranger-atac 2.1.0 count` × 9 samples; `--chemistry=ARC-v1`; ARC reference (`refdata-cellranger-arc-GRCh38-2020-A-2.0.0`) |
| `02_qc_filter.ipynb` | `snap.pp.filter_cells` — TSSe ≥ 6, `passed_filters + 1 ∈ [5000, 50000]`, `is__cell_barcode == 1` |
| `03_doublet_removal.ipynb` | 500-bp tile matrix, `select_features(n_features=250000)`, `snap.pp.scrublet`, `filter_doublets(probability_threshold=0.5)` |
| `04_clustering_annotation.ipynb` | `select_features(n_features=200000)`, Harmony (`batch="sample", max_iter_harmony=20`), kNN (k=10), Leiden (random_state=0), manual cell-type labels |

Note: `select_features` uses **250,000 features per sample** in step 03 and **200,000 features on the merged AnnDataSet** in step 04. See `s1_preprocessing/README.md`.

Input accessions are listed in Additional file 1, Table S1. Output: 44,262 nuclei across ten cell types.

## S2. Peak calling and construction of cell-type peak sets

| Script | Purpose |
|---|---|
| `01_make_tn5_bed.sh` | Fragment ends as 1-bp insertion events |
| `02_macs2_standard.sh` | MACS2, default q ≤ 0.05 |
| `03_macs2_permissive.sh` | MACS2, `-p 0.01` |
| `04_merge_standard_snapatac2.py` | `snap.tl.merge_peaks(half_width=250)` — 501-bp fixed-width per cell type (Standard) |
| `05_cluster_summits_3d.py` | 3D Euclidean clustering over (start, summit, end) with `--max_gap 70`, MM distance-weighted centroid (tol 1e-6, max 100 iters), asymmetric interval boundaries (Suggested) |

The `bedtools intersect -f 0.5 -r` deduplication for the Suggested workflow
moved to `s3_gkmqc/01_dedup_50pct.py` (its output is one of two gkmQC input
peak sets, together with the raw MM centroid).

## S3. Sequence-based quality assessment and trimming

**Two branches, on-disk fact — do not conflate.** Full comparison table in
`s3_gkmqc/README.md`.

| Script | Branch | Purpose |
|---|---|---|
| `01_dedup_50pct.py`            | A only        | 50 % reciprocal-overlap dedup on the raw MM centroid |
| `02_run_gkmqc_cpu.sh`          | A (Fig 4a)    | CPU `gkmqc.py evaluate -re 40` on the dedup'd set |
| `03_trim_first_below.py`       | A (Fig 4a)    | **first-below-0.7** rule, cap 40 subsets → 200,000 peaks; produces the Table S2 / Fig 4a counts |
| `04_run_gkmqc_cpu_for_ldsc.sh` | B (LDSC)      | CPU `gkmqc.py evaluate -re 350` on the raw MM centroid (manuscript run used GPU build; CPU here is equivalent, see `s3_gkmqc/README.md`) |
| `05_trim_max_above_for_ldsc.py`| B (LDSC)      | **max-above-0.7** rule; produces the hg38 BEDs feeding S5 |

gkmQC's null-sequence index is generated internally (UCSC hg38); no separate `build_null_index` script exists.

## S4. Comparison of peak sets

| Script | Purpose |
|---|---|
| `01_concordance.py` | Build ±150 bp summit windows (in-script) and compute asymmetric matched fractions via `bedtools intersect -a A -b B -u -sorted`; per-cell-type overlap TSV consumed by `figures/fig4_peak_overlap.py` |

## S5. Partitioned heritability

| Script | Purpose |
|---|---|
| `00_phase1_bed_ds.py` | Assemble hg38 per-cell-type BEDs for the Standard annotation (Suggested BEDs come from `s3_gkmqc/05_trim_max_above_for_ldsc.py`) |
| `01_munge_sumstats.sh` | `munge_sumstats.py` × 4 GWAS against the HapMap3 allele list |
| `02_liftover.sh` | **CrossMap 0.7.0** hg38 → hg19 (±1 kb extension is applied upstream in phase 1, in pandas — no `bedtools merge`) |
| `03_make_annot.sh` | `make_annot.py` per chromosome × cell type (single-annotation LD-score input, phase 4) |
| `03b_annot_overlap.sh` | `intersectBed -c` to append cell-type binary columns to baselineLD v2.2 (joint LD-score input, phase 7) |
| `04_ldscore.sh` | `ldsc.py --l2 --thin-annot --ld-wind-cm 1` (single-annot) |
| `04b_ldscore_overlap.sh` | `ldsc.py --l2 --ld-wind-cm 1` on the overlap annots (no `--thin-annot`) |
| `05_partitioned_h2.sh` | `ldsc.py --h2 --overlap-annot --print-coefficients` — main joint analysis (baselineLD v2.2 + 10 cell-type annots per workflow) |
| `06_h2cts_sensitivity.sh` | `ldsc.py --h2-cts` on all 20 (10 CT × 2 workflow) annotations |

GWAS sources are listed in Additional file 1, Table S3.

## Figures

Figure scripts live in `figures/`. `fig4_peak_overlap.py` builds Fig 4a from
the S4 concordance TSV. `fig3_ldsc_enrichment.py` and `fig3_tau_zscore.py`
build the heritability heatmaps from `LDSC_master_table.csv`, which is
produced by `figures/make_master_table.py` reading S5 outputs.

---

## Software

| Software | Version |
|---|---|
| Cell Ranger ATAC | 2.1.0 |
| SnapATAC2 | 2.8.0 |
| harmonypy | 0.0.10 |
| scanpy | 1.10.3 |
| MACS2 | 2.2.6 |
| bedtools | 2.31.1 |
| CrossMap | 0.7.0 |
| gkmQC | 1.0.0 |
| LDSC | 1.0.1 |

Four conda environments are required: Python 3.9.21 for SnapATAC2, a separate environment for MACS2, Python 3.7.12 for gkmQC, and Python 2.7.18 for LDSC. See `envs/`.

## Data availability

Raw snATAC-seq data are available from the ENCODE Project under the accessions in Additional file 1, Table S1. GWAS summary statistics are available from the publications in Additional file 1, Table S3. Processed peak sets are archived separately (see `results/README.md`).

## Citation

See `CITATION.cff`.

## License

BSD 3-Clause. See `LICENSE`.
