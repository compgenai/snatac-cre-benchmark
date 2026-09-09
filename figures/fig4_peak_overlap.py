#!/usr/bin/env python3
"""
Split v3 figure into two standalone files:
  - figure_panelA_only_..._stdSnapMerge_..._v3.{pdf,png}  (peak count bar chart)
  - figure_panelB_only_..._stdSnapMerge_..._v3.{pdf,png}  (AUROC curves, joined)

Same data and styling as
  make_figure_main_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered_v3.py
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Overrides:
#   FIG_OUT    — figures/figure_benchmark directory (source_data lives inside)
#   DATA_ROOT  — directory containing Heart_celltype_cell_counts.csv
OUT         = Path(os.environ.get("FIG_OUT",   "./figures/figure_benchmark"))
SRC         = OUT / "source_data"
COUNTS_CSV  = Path(os.environ.get("DATA_ROOT", "./output_so")) / "Heart_celltype_cell_counts.csv"

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

RANK_MAX    = 40
AUROC_FLOOR = 0.7
Y_LIM       = (0.5, 1.0)
COL_STD = "#7f7f7f"
COL_SUG = "#d62728"
COLOR_B = {"Standard": COL_STD, "Suggested dedup_50pct": COL_SUG}
METHOD_ORDER = ["Standard", "Suggested dedup_50pct"]

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "savefig.dpi":     300,
    "savefig.bbox":    "tight",
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})

# ---------- load data ----------
ov = pd.read_csv(SRC / "peak_overlap_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered.tsv", sep="\t")
ov_idx = ov.set_index("cell_type").loc[CELLTYPES]
auc_std = pd.read_csv(SRC / "auroc_curves_allcells_stdSnapMerge.tsv", sep="\t")
auc_std = auc_std[auc_std["rank"] <= RANK_MAX].copy()
auc_d50 = pd.read_csv(SRC / "auroc_curves_allcells_dedup50pct.tsv", sep="\t")
auc_d50 = auc_d50[auc_d50["rank"] <= RANK_MAX].copy()
auc = pd.concat([auc_std, auc_d50], ignore_index=True)
auc = auc[auc["method"].isin(METHOD_ORDER)].copy()
fb = (auc.dropna(subset=["first_below_0.7_rank"])
          .groupby(["method", "cell_type"])["first_below_0.7_rank"]
          .first().to_dict())


def fmt_k(n):
    if n >= 1000:
        return f"{n/1000:.0f}k"
    return f"{n}"


def draw_panel_A(ax):
    n_ct = len(CELLTYPES)
    y = np.arange(n_ct)[::-1].astype(float)
    bar_h = 0.36
    offset = bar_h / 2 + 0.02
    std_tot  = ov_idx["std_total"].values
    sug_tot  = ov_idx["sug_total"].values
    std_cpct = ov_idx["std_common_pct"].values
    sug_cpct = ov_idx["sug_common_pct"].values

    ax.barh(y + offset, std_tot, height=bar_h, color=COL_STD,
            edgecolor="black", lw=0.4, label="Standard snapMerge (all peaks)")
    ax.barh(y - offset, sug_tot, height=bar_h, color=COL_SUG,
            edgecolor="black", lw=0.4, label="Suggested dedup_50pct (AUC≥0.7, tie-aware)")

    xmax = max(std_tot.max(), sug_tot.max())
    pad  = xmax * 0.012
    for i in range(n_ct):
        ax.text(std_tot[i] + pad, y[i] + offset,
                f"n={std_tot[i]:,}  (common {std_cpct[i]:.1f}%)",
                va="center", ha="left", fontsize=6.8, color="black")
        ax.text(sug_tot[i] + pad, y[i] - offset,
                f"n={sug_tot[i]:,}  (common {sug_cpct[i]:.1f}%)",
                va="center", ha="left", fontsize=6.8, color="black")

    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[c] for c in CELLTYPES])
    for tick, ct in zip(ax.get_yticklabels(), CELLTYPES):
        if TIER[ct] == "rare":
            tick.set_color("crimson"); tick.set_fontweight("bold")

    ax.set_xlim(0, xmax * 1.32)
    xt = np.linspace(0, xmax, 6)
    xt = np.round(xt / 50000) * 50000
    ax.set_xticks(xt)
    ax.set_xticklabels([fmt_k(int(t)) for t in xt])
    ax.set_xlabel("Number of peaks")
    ax.set_title("Peak counts per cell type — Standard snapMerge vs Suggested dedup_50pct "
                 "(parentheses: % of own set that overlaps the other)")
    ax.legend(loc="lower right", frameon=False, ncol=1)
    ax.grid(axis="x", linestyle=":", alpha=0.4)


def draw_panel_B(axes_2x5):
    for i in range(2):
        for j in range(5):
            ct = CELLTYPES[i * 5 + j]
            ax = axes_2x5[i, j]
            sub_ct = auc[auc["cell_type"] == ct]
            ax.axhline(AUROC_FLOOR, color="black", lw=0.6, ls="--", alpha=0.8, zorder=1)
            for m in METHOD_ORDER:
                g = sub_ct[sub_ct["method"] == m].sort_values("rank")
                if g.empty: continue
                x = g["rank"].values; yv = g["auc"].values
                col = COLOR_B[m]
                cutoff = fb.get((m, ct), np.nan)
                if np.isnan(cutoff):
                    ax.plot(x, yv, "-", color=col, lw=1.6, label=m, zorder=3)
                else:
                    ms_solid = x < cutoff; ms_fade = x >= cutoff
                    if ms_solid.any():
                        xs, ys = x[ms_solid], yv[ms_solid]
                        if ms_fade.any():
                            xs = np.append(xs, x[ms_fade][0])
                            ys = np.append(ys, yv[ms_fade][0])
                        ax.plot(xs, ys, "-", color=col, lw=1.6, label=m, zorder=3)
                    if ms_fade.any():
                        xf, yf = x[ms_fade], yv[ms_fade]
                        ax.plot(xf, yf, "--", color=col, lw=1.2, alpha=0.45, zorder=3)
            title_color = "crimson" if TIER[ct] == "rare" else "black"
            ax.text(0.5, 0.97, LABEL[ct], transform=ax.transAxes,
                    ha="center", va="top", fontsize=8.5, color=title_color,
                    fontweight="bold" if TIER[ct] == "rare" else "normal")

    axes_2x5[0, 0].set_xlim(1, RANK_MAX)
    axes_2x5[0, 0].set_ylim(*Y_LIM)
    axes_2x5[0, 0].set_yticks(np.arange(0.5, 1.001, 0.1))
    for i in range(2):
        for j in range(5):
            ax = axes_2x5[i, j]
            if j != 0: ax.tick_params(labelleft=False)
            if i != 1: ax.tick_params(labelbottom=False)
    for j in range(1, 5):
        axes_2x5[1, j].tick_params(labelbottom=False)
    axes_2x5[0, 0].set_ylabel("AUC")
    axes_2x5[1, 0].set_xlabel("Rank bin (5,000 peaks each)")


# ---------- Panel A standalone ----------
figA, axA = plt.subplots(figsize=(13.0, 5.6))
draw_panel_A(axA)
pdfA = OUT / "figure_panelA_only_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered_v3.pdf"
pngA = OUT / "figure_panelA_only_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered_v3.png"
figA.savefig(pdfA); figA.savefig(pngA); plt.close(figA)
print(f"[saved] {pdfA}")
print(f"[saved] {pngA}")

# ---------- Panel B standalone (joined 2x5) ----------
figB = plt.figure(figsize=(13.0, 6.0))
gsB = GridSpec(2, 5, figure=figB, wspace=0, hspace=0)
axB = np.empty((2, 5), dtype=object)
for i in range(2):
    for j in range(5):
        if i == 0 and j == 0:
            ax = figB.add_subplot(gsB[i, j])
        else:
            ax = figB.add_subplot(gsB[i, j], sharex=axB[0, 0], sharey=axB[0, 0])
        axB[i, j] = ax
draw_panel_B(axB)

handles = [plt.Line2D([0], [0], color=COLOR_B[m], lw=2, label=m) for m in METHOD_ORDER]
figB.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
            bbox_to_anchor=(0.5, -0.02))
figB.suptitle("AUC curves — Standard snapMerge vs Suggested dedup_50pct "
              "(AUC≥0.7 tie-aware) — all cell types (ordered by cell count)",
              y=0.995, fontsize=11, fontweight="bold")

pdfB = OUT / "figure_panelB_only_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered_v3.pdf"
pngB = OUT / "figure_panelB_only_allcells_dedup50pct_auc07_tieaware_stdSnapMerge_fullname_ordered_v3.png"
figB.savefig(pdfB); figB.savefig(pngB); plt.close(figB)
print(f"[saved] {pdfB}")
print(f"[saved] {pngB}")
