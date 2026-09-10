# S2 · Peak calling and cell-type peak-set construction

Runs both workflows on the same per-(sample × cell type) pseudobulk Tn5
insertion tracks. Both start from per-pseudobulk MACS2 calls, but merge to
a per-cell-type peak set via different routes:

- **Standard workflow**: `macs2 callpeak` at the default threshold (q ≤ 0.05)
  per pseudobulk (outputs at `$DATA_ROOT/peakcalling_default/`), then
  `snap.tl.merge_peaks(half_width=250)` — 501-bp fixed-width iterative
  overlap removal per cell type.
- **Suggested workflow**: `macs2 callpeak -p 0.01` per pseudobulk (outputs at
  `$DATA_ROOT/peakcalling_p0.01/`), then 3D Euclidean summit clustering + MM
  distance-weighted centroid with asymmetric intervals to produce a
  per-cell-type peak set. The subsequent deduplication and AUC-trim that
  finalize the Fig 4a peak set live in **S3** — see `s3_gkmqc/README.md`.

## Pipeline layout

```
peakcalling_default/  ──→  04_merge_standard_snapatac2.py  ──→  Std_snapMerge_<CT>.narrowPeak   (Standard)
peakcalling_p0.01/    ──→  05_cluster_summits_3d.py        ──→  Improved_<CT>.narrowPeak         (Suggested)
```

## Scripts (run in order)

| # | File | Purpose |
|---|------|---------|
| 01 | `01_make_tn5_bed.sh`               | Extract 1-bp Tn5 insertion events from both fragment ends (no +4/-5 shift) |
| 02 | `02_macs2_standard.sh`             | MACS2, default q ≤ 0.05 (Standard) — per pseudobulk → `peakcalling_default/` |
| 03 | `03_macs2_permissive.sh`           | MACS2, `-p 0.01` (Suggested) — per pseudobulk → `peakcalling_p0.01/` |
| 04 | `04_merge_standard_snapatac2.py`   | `snap.tl.merge_peaks(half_width=250)` — 501-bp fixed-width per cell type (Standard) |
| 05 | `05_cluster_summits_3d.py`         | 3D Euclidean clustering over (start, summit, end) with `--max_gap 70`, MM distance-weighted centroid (tol 1e-6, max 100 iters), asymmetric interval boundaries (Suggested) |

## Inputs

- Per-(sample × cell type) fragment BEDs from S1 · 04 (`$DATA_ROOT/fragments/`)
- Tn5 insertion BEDs written by 01 (`$DATA_ROOT/tn5_insertions/`)

## Outputs (not tracked)

- Per-(sample × cell type) MACS2 narrowPeak:
  - Standard: `$DATA_ROOT/peakcalling_default/mcluster<CT>.<SAMPLE>/macs2_*/`
  - Suggested: `$DATA_ROOT/peakcalling_p0.01/mcluster<CT>.<SAMPLE>/macs2_*/`
- Per-cell-type narrowPeak (feeds S3):
  - Standard: `$OUT_ROOT/standard_snapMerge/Std_snapMerge_<CT>.narrowPeak` (501-bp fixed width)
  - Suggested: `$OUT_ROOT/gkmQC/MM/Improved_<CT>.narrowPeak` (raw MM centroid)

## Environments

- Steps 01, 04, 05: `envs/snapatac2.yml`
- Steps 02, 03: `envs/macs2.yml` (MACS2 2.2.6)
