#!/usr/bin/env python3
"""
3-way peak analysis: Standard-specific / Common / Suggested-specific.

Windows (same as figure v2, AUC>=0.7 tie-aware setup):
  Standard  = standard_snapMerge/Std_snapMerge_{ct}.narrowPeak, ALL peaks, ±150 bp around summit
  Suggested = gkmQC/Sug_dedup_50pct/Improved_{ct}.dedup_50pct.narrowPeak,
              AUC>=0.7 tie-aware cutoff (rank n_auc07_dedup50pct), ±150 bp around summit

Per cell type outputs (BED, 3 columns):
  {ct}_std_specific.bed  — Standard peaks with no overlap in Suggested
  {ct}_common_std.bed    — Standard peaks that overlap Suggested (Std-side view of common)
  {ct}_common_sug.bed    — Suggested peaks that overlap Standard (Sug-side view of common)
  {ct}_sug_specific.bed  — Suggested peaks with no overlap in Standard

Aggregate output:
  3way_counts.tsv  — per-ct counts + % for all 3 categories, both views
  3way_figure.pdf/.png  — stacked bar visualization
"""
import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------- config ----------
# Overrides:
#   DATA_ROOT — directory containing Heart_celltype_cell_counts.csv (S1 output)
#   OUT_ROOT  — directory containing gkmQC/ and standard_snapMerge/
#   FIG_OUT   — output directory for windowed BEDs, overlap TSVs, figures
DATA_ROOT   = Path(os.environ.get("DATA_ROOT", "./output_so"))
OUT_ROOT    = Path(os.environ.get("OUT_ROOT",  "./output_seung"))
FIG_OUT     = Path(os.environ.get("FIG_OUT",   "./figures/figure_benchmark"))

SEUNG       = OUT_ROOT / "gkmQC"
SNAPMERGE   = OUT_ROOT / "standard_snapMerge"
OUT         = FIG_OUT
TMP_WIN     = OUT / "tmp_overlap_dedup50pct_auc07_tieaware_stdSnapMerge"
WORK        = OUT / "3way_stdSnapMerge_dedup50pct"
BEDS        = WORK / "beds"
SRC         = OUT / "source_data"
COUNTS_CSV  = DATA_ROOT / "Heart_celltype_cell_counts.csv"
for d in (TMP_WIN, WORK, BEDS):
    d.mkdir(parents=True, exist_ok=True)

BEDTOOLS = os.environ.get("BEDTOOLS", "bedtools")

WIN = 150
COL_STD_SPEC = "#7f7f7f"
COL_COMMON   = "#2ca02c"
COL_SUG_SPEC = "#d62728"

# Cell-type order: total cell count descending
_counts = pd.read_csv(COUNTS_CSV)
_counts = _counts[_counts["Cell_Type"] != "Total"].copy()
_counts["Total_Cells"] = _counts["Total_Cells"].astype(int)
_counts = _counts.sort_values("Total_Cells", ascending=False).reset_index(drop=True)
CELLTYPES = _counts["Cell_Type"].tolist()

LABEL = {ct: ct.replace("_", " ") for ct in CELLTYPES}
TIER = {
    "Ventricular_Cardiomyocytes": "common",
    "Primitive_Endoderm":         "other",
    "Fibroblasts":                "common",
    "Atrial_Cardiomyocytes":      "other",
    "Endothelial":                "other",
    "Smooth_Muscle":              "other",
    "Macrophages":                "other",
    "Nervous_Cells":              "rare",
    "Myofibroblasts":             "other",
    "Trophectoderm":              "rare",
}


# ---------- helpers ----------
def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout

def wc(p):
    p = Path(p)
    if not p.exists() or p.stat().st_size == 0:
        return 0
    return int(subprocess.check_output(["wc", "-l", str(p)]).split()[0])


def build_std_window(ct):
    """Standard (snapMerge) ±WIN bp around summit; identical to v2 script cache."""
    out = TMP_WIN / f"{ct}_std_snapMerge_{2*WIN}bp.bed"
    if out.exists() and out.stat().st_size > 0:
        return out
    df = pd.read_csv(SNAPMERGE / f"Std_snapMerge_{ct}.narrowPeak", sep="\t", header=None)
    summit = df[1] + df[9]
    bed = pd.DataFrame({0: df[0],
                        "s": (summit - WIN).clip(lower=0),
                        "e": summit + WIN}).sort_values([0, "s"])
    bed.to_csv(out, sep="\t", header=False, index=False)
    return out


def build_sug_window(ct, n_sug):
    """Suggested dedup_50pct, AUC>=0.7 tie-aware; identical to v2 script cache."""
    out = TMP_WIN / f"{ct}_sug_tieaware_{2*WIN}bp.bed"
    if out.exists() and out.stat().st_size > 0:
        return out
    df = pd.read_csv(SEUNG / "Sug_dedup_50pct" / f"Improved_{ct}.dedup_50pct.narrowPeak",
                     sep="\t", header=None)
    df_sorted = df.sort_values(7, ascending=False).reset_index(drop=True)
    cutoff = float(df_sorted.iloc[n_sug - 1, 7])
    df_filt = df_sorted[df_sorted[7] >= cutoff].copy()
    summit = df_filt[1] + df_filt[9]
    df_filt["s"] = (summit - WIN).clip(lower=0)
    df_filt["e"] = summit + WIN
    bed = df_filt[[0, "s", "e"]].sort_values([0, "s"])
    bed.to_csv(out, sep="\t", header=False, index=False)
    return out


# ---------- load AUC>=0.7 targets ----------
red = pd.read_csv(SRC / "redundancy_with_auroc_filter_50pct.tsv",
                  sep="\t").set_index("cell_type")


# ---------- 3-way BED + counts ----------
print("== Building 3-way BED files per cell type ==")
rows = []
for ct in CELLTYPES:
    n_sug = int(red.loc[ct, "n_auc07_dedup50pct"])

    std_w = build_std_window(ct)
    sug_w = build_sug_window(ct, n_sug)

    std_spec = BEDS / f"{ct}_std_specific.bed"
    com_std  = BEDS / f"{ct}_common_std.bed"
    com_sug  = BEDS / f"{ct}_common_sug.bed"
    sug_spec = BEDS / f"{ct}_sug_specific.bed"

    # Standard-specific: std ∖ sug
    sh(f"{BEDTOOLS} intersect -a {std_w} -b {sug_w} -v -sorted > {std_spec}")
    # Standard-side common: std ∩ sug (std peaks with any overlap in sug)
    sh(f"{BEDTOOLS} intersect -a {std_w} -b {sug_w} -u -sorted > {com_std}")
    # Suggested-side common: sug ∩ std
    sh(f"{BEDTOOLS} intersect -a {sug_w} -b {std_w} -u -sorted > {com_sug}")
    # Suggested-specific: sug ∖ std
    sh(f"{BEDTOOLS} intersect -a {sug_w} -b {std_w} -v -sorted > {sug_spec}")

    n_std_total = wc(std_w)
    n_sug_total = wc(sug_w)
    n_std_spec  = wc(std_spec)
    n_com_std   = wc(com_std)
    n_com_sug   = wc(com_sug)
    n_sug_spec  = wc(sug_spec)

    assert n_std_spec + n_com_std == n_std_total, f"{ct}: std split mismatch"
    assert n_sug_spec + n_com_sug == n_sug_total, f"{ct}: sug split mismatch"

    rows.append({
        "cell_type":         ct,
        "label":             LABEL[ct],
        "tier":              TIER[ct],
        "n_sug_target":      n_sug,
        "std_total":         n_std_total,
        "std_specific":      n_std_spec,
        "common_std_side":   n_com_std,
        "common_sug_side":   n_com_sug,
        "sug_specific":      n_sug_spec,
        "sug_total":         n_sug_total,
        "std_specific_pct":  round(100*n_std_spec/n_std_total, 2) if n_std_total else 0,
        "common_std_pct":    round(100*n_com_std/n_std_total, 2)  if n_std_total else 0,
        "common_sug_pct":    round(100*n_com_sug/n_sug_total, 2)  if n_sug_total else 0,
        "sug_specific_pct":  round(100*n_sug_spec/n_sug_total, 2) if n_sug_total else 0,
    })
    print(f"  {ct:30s} std_total={n_std_total:>7,}  "
          f"std_spec={n_std_spec:>6,}  com(std/sug)={n_com_std:>6,}/{n_com_sug:>6,}  "
          f"sug_spec={n_sug_spec:>6,}  sug_total={n_sug_total:>7,}")

df = pd.DataFrame(rows)
tsv = WORK / "3way_counts.tsv"
df.to_csv(tsv, sep="\t", index=False)
print(f"[saved] {tsv}")


# ---------- Figure ----------
plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "savefig.dpi":     300,
    "savefig.bbox":    "tight",
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})

df_idx = df.set_index("cell_type").loc[CELLTYPES]
n_ct = len(CELLTYPES)
y = np.arange(n_ct)[::-1].astype(float)
bar_h = 0.36
offset = bar_h / 2 + 0.02

fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                         gridspec_kw={"width_ratios": [1, 1], "wspace": 0.35})

# --- Panel A: absolute counts (paired bars, Standard top / Suggested bottom) ---
ax = axes[0]
for i, ct in enumerate(CELLTYPES):
    r = df_idx.loc[ct]
    # Standard bar (top of pair)
    ax.barh(y[i] + offset, r["common_std_side"],       height=bar_h,
            color=COL_COMMON, edgecolor="black", lw=0.4)
    ax.barh(y[i] + offset, r["std_specific"], left=r["common_std_side"],
            height=bar_h, color=COL_STD_SPEC, edgecolor="black", lw=0.4)
    # Suggested bar (bottom of pair)
    ax.barh(y[i] - offset, r["common_sug_side"],       height=bar_h,
            color=COL_COMMON, edgecolor="black", lw=0.4)
    ax.barh(y[i] - offset, r["sug_specific"], left=r["common_sug_side"],
            height=bar_h, color=COL_SUG_SPEC, edgecolor="black", lw=0.4)

xmax = max(df_idx[["std_total", "sug_total"]].max())
pad  = xmax * 0.012
for i, ct in enumerate(CELLTYPES):
    r = df_idx.loc[ct]
    ax.text(r["std_total"] + pad, y[i] + offset,
            f"Std n={r['std_total']:,} (spec {r['std_specific_pct']:.1f}%)",
            va="center", ha="left", fontsize=6.5)
    ax.text(r["sug_total"] + pad, y[i] - offset,
            f"Sug n={r['sug_total']:,} (spec {r['sug_specific_pct']:.1f}%)",
            va="center", ha="left", fontsize=6.5)

ax.set_yticks(y)
ax.set_yticklabels([LABEL[c] for c in CELLTYPES])
for tick, ct in zip(ax.get_yticklabels(), CELLTYPES):
    if TIER[ct] == "rare":
        tick.set_color("crimson"); tick.set_fontweight("bold")

ax.set_xlim(0, xmax * 1.38)
ax.set_xlabel("Number of peaks")
ax.set_title("Peak counts — common vs specific  (top: Standard snapMerge, bottom: Suggested)")
ax.grid(axis="x", linestyle=":", alpha=0.4)

# --- Panel B: normalized (each row = 100%, per-side stacked) ---
ax = axes[1]
for i, ct in enumerate(CELLTYPES):
    r = df_idx.loc[ct]
    std_com_p = 100 * r["common_std_side"] / r["std_total"]
    std_spec_p = 100 * r["std_specific"]   / r["std_total"]
    sug_com_p = 100 * r["common_sug_side"] / r["sug_total"]
    sug_spec_p = 100 * r["sug_specific"]   / r["sug_total"]
    ax.barh(y[i] + offset, std_com_p,  height=bar_h, color=COL_COMMON, edgecolor="black", lw=0.4)
    ax.barh(y[i] + offset, std_spec_p, left=std_com_p, height=bar_h, color=COL_STD_SPEC, edgecolor="black", lw=0.4)
    ax.barh(y[i] - offset, sug_com_p,  height=bar_h, color=COL_COMMON, edgecolor="black", lw=0.4)
    ax.barh(y[i] - offset, sug_spec_p, left=sug_com_p, height=bar_h, color=COL_SUG_SPEC, edgecolor="black", lw=0.4)
    # % labels centered inside common segment
    if std_com_p >= 8:
        ax.text(std_com_p/2, y[i] + offset, f"{std_com_p:.0f}%",
                va="center", ha="center", fontsize=6.5, color="white", fontweight="bold")
    if sug_com_p >= 8:
        ax.text(sug_com_p/2, y[i] - offset, f"{sug_com_p:.0f}%",
                va="center", ha="center", fontsize=6.5, color="white", fontweight="bold")

ax.set_yticks(y)
ax.set_yticklabels([LABEL[c] for c in CELLTYPES])
for tick, ct in zip(ax.get_yticklabels(), CELLTYPES):
    if TIER[ct] == "rare":
        tick.set_color("crimson"); tick.set_fontweight("bold")
ax.set_xlim(0, 100)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
ax.set_xlabel("Fraction within each set")
ax.set_title("Composition of each set  (top: Standard, bottom: Suggested)")
ax.grid(axis="x", linestyle=":", alpha=0.4)

# --- shared legend ---
handles = [
    Patch(facecolor=COL_STD_SPEC, edgecolor="black", label="Standard-specific"),
    Patch(facecolor=COL_COMMON,   edgecolor="black", label="Common (overlap)"),
    Patch(facecolor=COL_SUG_SPEC, edgecolor="black", label="Suggested-specific"),
]
fig.legend(handles=handles, loc="lower center", ncol=3,
           frameon=False, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("3-way peak overlap — Standard snapMerge vs Suggested dedup_50pct (AUC≥0.7 tie-aware)",
             y=1.02, fontsize=12, fontweight="bold")

pdf = WORK / "3way_figure.pdf"
png = WORK / "3way_figure.png"
plt.savefig(pdf); plt.savefig(png); plt.close()
print(f"[saved] {pdf}")
print(f"[saved] {png}")


# ---------- summary ----------
print("\n== 3-way counts ==")
show_cols = ["cell_type", "tier", "std_total", "std_specific", "common_std_side",
             "common_sug_side", "sug_specific", "sug_total"]
print(df[show_cols].to_string(index=False))
