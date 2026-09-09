# S4 · Peak-set comparison (Fig 4a source)

Computes asymmetric matched fractions (concordance %) between the Standard
and Suggested cell-type peak sets on **±150 bp summit windows** (300-bp
intervals). Feeds the per-cell-type overlap TSV that the Fig 4a bar chart
in `figures/fig4_peak_overlap.py` consumes.

## Script

| # | File | Purpose |
|---|------|---------|
| 01 | `01_concordance.py` | Build ±150 bp summit windows for both sets, count peaks in each set that overlap at least one peak of the other via `bedtools intersect -a A -b B -u -sorted`, emit per-cell-type TSV |

The ±150 bp summit-window construction is embedded in `01_concordance.py`
(`build_std_window` / `build_sug_window` helpers, `WIN = 150`). There is no
separate `01_summit_windows.sh` — the two operations are one script.

The 2-kb-flanked intervals mentioned in some earlier working notes are for
**LDSC annotation** (S5), not for Fig 4a. Do not conflate the two.

## Inputs

- Standard: `output_seung/standard_snapMerge/Std_snapMerge_<CT>.narrowPeak`
- Suggested: `output_seung/gkmQC/Sug_dedup_50pct/Improved_<CT>.dedup_50pct.narrowPeak`
  (Branch A peak set, signalValue-tie-aware filtered to
  `n_auc07_dedup50pct` peaks from `s3_gkmqc/03_trim_first_below.py`)

## Outputs (not tracked)

- Per-cell-type overlap TSV
  (`peak_overlap_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_*.tsv`),
  consumed by `figures/fig4_peak_overlap.py`.

## Environment

`envs/snapatac2.yml` (bedtools 2.31.1, pandas 2.0.3).
