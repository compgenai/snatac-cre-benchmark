"""Enrichment heatmap by phenotype: 4 phenotypes x (DS | MM).  [v3]

Cell-type row order is by total cell count (descending) from
Heart_celltype_cell_counts.csv (intersected with cell types present in
LDSC_master_table.csv).
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# Overrides:
#   LDSC_WORK  — directory containing LDSC_master_table.csv
#   DATA_ROOT  — directory containing Heart_celltype_cell_counts.csv
ROOT       = Path(os.environ.get("LDSC_WORK", "./LDSC_work/v3"))
OUT        = ROOT / "figures"
COUNTS_CSV = Path(os.environ.get("DATA_ROOT", "./output_so")) / "Heart_celltype_cell_counts.csv"
OUT.mkdir(exist_ok=True)
df = pd.read_csv(ROOT / "LDSC_master_table.csv")

PHENO_ORDER = [
    "Aragam_2022_CAD_primary",
    "Nauffal_2022_QT",
    "Ntalla_2020_PR_EUR",
    "Roselli_2025_AF_common",
]
PHENO_LABEL = {
    "Aragam_2022_CAD_primary": "CAD\n(Aragam 2022)",
    "Nauffal_2022_QT":         "QT interval\n(Nauffal 2022)",
    "Ntalla_2020_PR_EUR":      "PR interval\n(Ntalla 2020)",
    "Roselli_2025_AF_common":  "AF\n(Roselli 2025)",
}
METHODS = ["DS", "MM"]
METHOD_LABEL = {"DS": "Standard", "MM": "Suggested"}

_counts = pd.read_csv(COUNTS_CSV)
_counts = _counts[_counts["Cell_Type"] != "Total"].copy()
_counts["Total_Cells"] = _counts["Total_Cells"].astype(int)
_counts = _counts.sort_values("Total_Cells", ascending=False).reset_index(drop=True)
_present = set(df["Celltype"].unique())
ct_order = [ct for ct in _counts["Cell_Type"].tolist() if ct in _present]
print("Cell-type order (by cell count desc, intersected with master table):")
for ct in ct_order:
    n = int(_counts.loc[_counts["Cell_Type"] == ct, "Total_Cells"].iloc[0])
    print(f"  {ct:30s} n={n:,}")

def stars(p):
    if pd.isna(p): return ""
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return ""

vmin = df["Enrichment"].quantile(0.02)
vmax = df["Enrichment"].quantile(0.98)
vmax = max(vmax, 1.0)
vmin = min(vmin, 0.0)

fig, axes = plt.subplots(1, 4, figsize=(13.5, 5.2),
                         gridspec_kw={"wspace": 0.45})

cmap = plt.cm.RdBu_r
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

for j, pheno in enumerate(PHENO_ORDER):
    ax = axes[j]
    sub = df[df.Phenotype == pheno]
    mat = np.full((len(ct_order), 2), np.nan)
    pmat = np.full((len(ct_order), 2), np.nan)
    for i, ct in enumerate(ct_order):
        for k, mthd in enumerate(METHODS):
            row = sub[(sub.Celltype == ct) & (sub.Category == mthd)]
            if len(row):
                mat[i, k]  = row["Enrichment"].iloc[0]
                pmat[i, k] = row["Enrichment_p"].iloc[0]

    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)
    for i in range(len(ct_order)):
        for k in range(2):
            if np.isnan(mat[i, k]):
                continue
            txt = f"{mat[i,k]:.1f}{stars(pmat[i,k])}"
            rgba = cmap(norm(mat[i, k]))
            brightness = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            tc = "white" if brightness < 0.5 else "black"
            ax.text(k, i, txt, ha="center", va="center",
                    fontsize=7.5, color=tc)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([METHOD_LABEL[m] for m in METHODS],
                       fontsize=9, rotation=0)
    if j == 0:
        ax.set_yticks(range(len(ct_order)))
        ax.set_yticklabels([ct.replace("_", " ") for ct in ct_order],
                           fontsize=9)
    else:
        ax.set_yticks([])
    ax.set_title(PHENO_LABEL[pheno], fontsize=10, pad=8)
    for x in [0.5]:
        ax.axvline(x, color="white", lw=1.2)
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ct_order), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.4)
    ax.tick_params(which="minor", length=0)

cax = fig.add_axes([0.92, 0.18, 0.014, 0.66])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
cb.set_label("Enrichment", fontsize=10)
cb.ax.tick_params(labelsize=8)

fig.suptitle("LDSC Partitioned h$^2$ Enrichment - Standard vs Suggested (v3, ordered by cell count)",
             fontsize=12, y=1.00)
fig.text(0.5, -0.01,
         "*** p<0.001   ** p<0.01   * p<0.05",
         ha="center", fontsize=8.5, style="italic")

plt.savefig(OUT / "fig_enrichment_heatmap_by_phenotype_ordered_by_count.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(OUT / "fig_enrichment_heatmap_by_phenotype_ordered_by_count.png",
            bbox_inches="tight", dpi=200)
print("saved:", OUT / "fig_enrichment_heatmap_by_phenotype_ordered_by_count.pdf")
