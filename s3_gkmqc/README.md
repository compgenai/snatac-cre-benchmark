# S3 · Sequence-based quality assessment (gkmQC) and trimming

**Two branches coexist**, on two different peak sets, with two different
AUC-trim rules. Both are on-disk fact — neither is derived from the other.
Downstream consumers differ: Branch A produces the peak-count series shown
in Table S2 / Fig 4a; Branch B produces the input to the LDSC Suggested
annotation.

The repository standardizes both branches on the **CPU** gkmQC build
(`/data/pipelines/gkmqc/bin/gkmqc.py`). The LDSC-branch run that generated
the manuscript's numbers originally used the GPU build (`gkmqc-gpu.py
-P slurm_gpu`) for wall-clock reasons; the two builds produce equivalent
`.eval.out` AUC profiles, so switching to CPU affects wall time only.

## Branch comparison

| Property | Branch A (Fig 4a / Table S2) | Branch B (LDSC Suggested annotation) |
|---|---|---|
| Input peak set | **After** 50 % reciprocal dedup (`Improved_<CT>.dedup_50pct.narrowPeak`) | **Before** dedup — raw MM 3D-centroid output (`Improved_<CT>.narrowPeak`) |
| Dedup applied? | Yes — via `01_dedup_50pct.py`                                                | No                                                                          |
| gkmQC build    | CPU `gkmqc.py`                                                               | CPU `gkmqc.py` (was `gkmqc-gpu.py` in manuscript run)                       |
| gkmQC wrapper  | `02_run_gkmqc_cpu.sh`                                                        | `04_run_gkmqc_cpu_for_ldsc.sh`                                              |
| `-re` (max subsets) | **40**                                                                  | **350**                                                                     |
| `-@` (threads) | 100                                                                          | 100                                                                         |
| AUC-trim rule  | **first-below 0.7**, cap 40 subsets → 200,000 peaks                          | **max-above 0.7** (largest subset index whose mean AUROC still exceeds 0.7) |
| Trim script    | `03_trim_first_below.py`                                                     | `05_trim_max_above_for_ldsc.py`                                             |
| Trim key line  | `usable_rank = int(first_below_0.7_rank) - 1; final_peak = usable_rank * 5000` | `last_topN = int(ed.loc[ed.auc > 0.7, "topN"].max())`                       |
| Non-monotone effect | 3 cell types (ACM, Endothelial, PE) end with 15,000–25,000 fewer peaks under this rule than under max-above | Retains peaks past a temporary dip if AUROC recovers                        |
| Output         | Peak-count target column in `redundancy_with_auroc_filter_50pct.tsv` (values shown in Table S2) | Per-cell-type hg38 BED at `${LDSC_WORK}/phase1_bed_hg38/MM_<CT>.bed` (feeds S5) |

For a full per-cell-type crossings table and Δ-peak analysis, see the
manuscript's Supplementary Methods audit
(`SUPP_CONFIRM_ROUND3.md` §Q1 / `SUPP_CONFIRM_ROUND4.md` §3).

## null-sequence index

**No user-authored `01_build_null_index.sh` exists.** gkmQC generates its
own GC- and repeat-matched negative sequences internally via
`/data/pipelines/gkmqc/scripts/seqs_nullgen.py` (invoked as part of
`gkmqc.py evaluate`), sampling from the UCSC hg38 assembly. Nothing to run
separately.

## Branch A (Fig 4a / Table S2) — run order

1. `01_dedup_50pct.py <input>.narrowPeak <output>.narrowPeak <tmpdir>` — 50 %
   reciprocal-overlap dedup on the raw MM centroid narrowPeak.
2. `02_run_gkmqc_cpu.sh <CELLTYPE>` — CPU gkmQC on the dedup'd set (`-re 40`).
3. `03_trim_first_below.py` — reads `.eval.out`, applies first-below-0.7
   rule, writes `n_auc07_dedup50pct` counts into the redundancy summary
   consumed by the Fig 4a builder in `s4_peak_comparison/`.

## Branch B (LDSC annotation) — run order

1. `04_run_gkmqc_cpu_for_ldsc.sh <CELLTYPE>` — CPU gkmQC on the raw MM
   centroid narrowPeak (`-re 350`).
2. `05_trim_max_above_for_ldsc.py` — reads `.eval.out`, applies max-above-0.7
   rule, filters by MACS2 -log10 pValue cutoff, extends ±1 kb, writes
   `${LDSC_WORK}/phase1_bed_hg38/MM_<CT>.bed` for S5.

## Env-var overrides

Every script in this folder reads its input/output roots and tool paths
from environment variables with sensible defaults. Common overrides:

```
export OUT_ROOT=/your/path/output_seung          # gkmQC/MM, gkmQC/Sug_dedup_50pct live here
export LDSC_WORK=/your/path/LDSC/v3              # Branch B trim script writes phase1_bed_hg38/
export GKMQC_ENV=/your/gkmqc/conda-env
export GKMQC_BIN=/your/gkmqc/bin/gkmqc.py
export BEDTOOLS=/your/bin/bedtools               # dedup step uses this
```

## Environment

`envs/gkmqc.yml` (Python 3.7.12, numpy 1.21.6, scikit-learn 1.0.2,
scipy 1.7.3, pyfasta 0.5.2, gkmQC v1.0.0). The dedup step (`01_dedup_50pct.py`)
also invokes `bedtools 2.31.1` from `envs/snapatac2.yml`.
