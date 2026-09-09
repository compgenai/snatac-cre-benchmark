#!/usr/bin/env python3
"""
Step 1: Add AUROC>=0.7 cutoff peak counts to redundancy summary.

For each (method, cell_type) get first_below_0.7_rank from auroc_summary_per_celltype.tsv,
compute usable_rank = (first_below_0.7_rank - 1), and final_peak = usable_rank * 5000.

If first_below_0.7_rank is NaN (= AUROC stays >=0.7 up to rank<=40), use 40 (analysis cap).

Outputs:
  analysis_redundancy/redundancy_with_auroc_filter.tsv
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np

# Override with:  export OUT_ROOT=/path/to/output_seung
BASE = Path(os.environ.get("OUT_ROOT", "./output_seung")) / "gkmQC"

red = pd.read_csv(BASE / "analysis_redundancy" / "redundancy_summary.tsv", sep="\t")
auc = pd.read_csv(BASE / "analysis_AUROC" / "auroc_summary_per_celltype.tsv", sep="\t")

fb = (auc.dropna(subset=["first_below_0.7_rank"])
         .groupby(["method", "cell_type"])["first_below_0.7_rank"]
         .first().to_dict())

RANK_MAX = 40
PEAKS_PER_RANK = 5000

def usable_rank(method, ct):
    val = fb.get((method, ct), np.nan)
    if pd.isna(val):
        return RANK_MAX
    return int(val) - 1

def filtered_peak_count(method, ct):
    return usable_rank(method, ct) * PEAKS_PER_RANK

HIGHLIGHT = {"Ventricular_Cardiomyocytes": "common",
             "Fibroblasts":                "common",
             "Nervous_Cells":              "rare",
             "Trophectoderm":              "rare"}

out = []
for _, row in red.iterrows():
    ct = row["cell_type"]
    out.append({
        "cell_type": ct,
        "tier": HIGHLIGHT.get(ct, "other"),
        # raw counts
        "standard_peaks":        int(row["standard_peaks"]),
        "suggested_raw_peaks":   int(row["suggested_raw_peaks"]),
        "suggested_dedup600":    int(row["dedup_600bp_peaks"]),
        # first_below_0.7 rank
        "fb_standard":           fb.get(("Standard", ct), np.nan),
        "fb_sug_raw":            fb.get(("Suggested raw", ct), np.nan),
        "fb_sug_dedup600":       fb.get(("Suggested dedup_600bp", ct), np.nan),
        # usable rank (last rank with AUROC >= 0.7, capped at 40)
        "usable_rank_standard":     usable_rank("Standard", ct),
        "usable_rank_sug_raw":      usable_rank("Suggested raw", ct),
        "usable_rank_sug_dedup600": usable_rank("Suggested dedup_600bp", ct),
        # final peak count after AUROC>=0.7 filter (= usable_rank * 5000)
        "auroc07_standard":     filtered_peak_count("Standard", ct),
        "auroc07_sug_raw":      filtered_peak_count("Suggested raw", ct),
        "auroc07_sug_dedup600": filtered_peak_count("Suggested dedup_600bp", ct),
    })

df = pd.DataFrame(out)
out_path = BASE / "analysis_redundancy" / "redundancy_with_auroc_filter.tsv"
df.to_csv(out_path, sep="\t", index=False)
print(f"[saved] {out_path}\n")

# Compact table sorted by standard peak count desc
df = df.sort_values("standard_peaks", ascending=False)

cols_compact = [
    "cell_type", "tier",
    "standard_peaks", "auroc07_standard",
    "suggested_raw_peaks", "auroc07_sug_raw",
    "suggested_dedup600", "auroc07_sug_dedup600",
]
print("== Peak counts: BEFORE vs AFTER AUROC>=0.7 filter ==")
print(df[cols_compact].to_string(index=False))

print("\n== Usable rank (last rank with AUROC>=0.7, capped at 40) ==")
print(df[["cell_type", "tier", "fb_standard", "fb_sug_raw", "fb_sug_dedup600",
         "usable_rank_standard", "usable_rank_sug_raw", "usable_rank_sug_dedup600"]]
      .to_string(index=False))

# Highlight 4-cell view
print("\n== Highlight 4 cell types (re-ordered) ==")
hl_order = ["Ventricular_Cardiomyocytes", "Fibroblasts", "Nervous_Cells", "Trophectoderm"]
df_hl = df.set_index("cell_type").loc[hl_order].reset_index()
print(df_hl[cols_compact].to_string(index=False))
