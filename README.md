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
| `01_cellranger_atac.sh` | Alignment to GRCh38, ARC-v1 chemistry |
| `02_qc_filter.py` | TSSe ≥ 6, 5,000–50,000 unique fragments, uniform thresholds across samples |
| `03_doublet_removal.py` | `snap.pp.scrublet` + `snap.pp.filter_doublets` |
| `04_clustering_annotation.py` | 500-bp tile matrix, Harmony, Leiden, manual cell-type assignment |

Input accessions are listed in Additional file 1, Table S1. Output: 44,262 nuclei across ten cell types.

## S2. Peak calling and construction of cell-type peak sets

| Script | Purpose |
|---|---|
| `01_make_tn5_bed.sh` | Fragment ends as 1-bp insertion events |
| `02_macs2_standard.sh` | MACS2, q ≤ 0.05 |
| `03_macs2_permissive.sh` | MACS2, p ≤ 0.01 |
| `04_merge_standard_snapatac2.py` | `snap.tl.merge_peaks(half_width=250)` |
| `05_cluster_summits_3d.py` | Euclidean clustering over (start, summit, end), 70 bp |
| `06_weiszfeld_centroid.py` | L1 geometric median, tol 1e-6, max 100 iterations |
| `07_dedup_reciprocal50.sh` | `bedtools intersect -f 0.5 -r`, lower-ranked member discarded |

Both workflows receive the same per-sample input.

## S3. Sequence-based quality assessment and trimming

| Script | Purpose |
|---|---|
| `01_build_null_index.sh` | Null-sequence index built locally from UCSC hg38 |
| `02_run_gkmqc.sh` | gkmQC `evaluate`, subsets of 5,000 peaks, 600-bp windows |
| `03_auc_trim.py` | Trim at the **first subset falling below** mean AUROC 0.7, tie-aware |

> The trimming rule is first-below-0.7. Record here which gkmQC build was used
> (CPU or GPU) and whether the input was the pre- or post-deduplication peak set,
> since the two are not interchangeable.

## S4. Comparison of peak sets

| Script | Purpose |
|---|---|
| `01_summit_windows.sh` | Summit-centered ±150 bp windows |
| `02_concordance.py` | Asymmetric matched fractions via `bedtools intersect -u` |

## S5. Partitioned heritability

| Script | Purpose |
|---|---|
| `01_munge_sumstats.sh` | `munge_sumstats.py` against the HapMap3 allele list |
| `02_liftover_extend.sh` | CrossMap hg38 → hg19, ±1 kb extension |
| `03_make_annot.sh` | Annotations appended to baselineLD v2.2 (1000G EUR Phase 3) |
| `04_ldscore.sh` | LD scores per chromosome, 1 cM window |
| `05_partitioned_h2.sh` | LDSC with `--overlap-annot --print-coefficients` |

GWAS sources are listed in Additional file 1, Table S3.

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
