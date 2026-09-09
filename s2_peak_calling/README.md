# S2 · Peak calling and cell-type peak-set construction

Runs both workflows (Standard MACS2 q ≤ 0.05 and Suggested MACS2 p ≤ 0.01)
on the same per-sample × per-cell-type pseudobulk Tn5 insertion tracks, then
constructs the corresponding cell-type peak sets:

- **Standard workflow**: `snap.tl.merge_peaks(half_width=250)` → 501-bp
  fixed-width peaks per cell type.
- **Suggested workflow**: 3D Euclidean summit clustering + MM
  distance-weighted centroid with asymmetric intervals (per cell type). The
  deduplication and AUC-trim that turn this set into the final Fig 4a peak
  set live in **S3** — see `s3_gkmqc/README.md`.

## Scripts (run in order)

| # | File | Purpose |
|---|------|---------|
| 01 | `01_make_tn5_bed.sh`               | Extract 1-bp Tn5 insertion events from both fragment ends (no +4/-5 shift) |
| 02 | `02_macs2_standard.sh`             | MACS2, default q ≤ 0.05 (Standard) |
| 03 | `03_macs2_permissive.sh`           | MACS2, `-p 0.01` (Suggested) |
| 04 | `04_merge_standard_snapatac2.py`   | `snap.tl.merge_peaks(half_width=250)` on Standard per-sample narrowPeak |
| 05 | `05_cluster_summits_3d.py`         | 3D Euclidean clustering (start, summit, end) with `--max_gap 70`, MM distance-weighted centroids (tol 1e-6, max 100 iters), asymmetric interval boundaries; single script fills both README slots `05_cluster_summits_3d.py` and `06_weiszfeld_centroid.py` |

## Inputs

- Post-annotation `AnnDataSet` from S1
- Per-cell-type / per-sample fragment BEDs (extracted by 01)

## Outputs (not tracked)

- Per-(sample × cell type) MACS2 narrowPeak (both q and p variants)
- Per-cell-type merged narrowPeak:
  - Standard: `Std_snapMerge_<CT>.narrowPeak` (501-bp fixed width)
  - Suggested: `Improved_<CT>.narrowPeak` (raw MM centroid — feeds S3 dedup)

## Environments

- Steps 01, 04, 05: `envs/snapatac2.yml`
- Steps 02, 03: `envs/macs2.yml` (MACS2 2.2.6)
