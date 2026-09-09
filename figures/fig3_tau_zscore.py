"""Coefficient z-score heatmap by phenotype.  [v3]"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# Override:  export LDSC_WORK=/path/to/LDSC/v3
ROOT = Path(os.environ.get("LDSC_WORK", "./LDSC_work/v3"))
OUT  = ROOT / "figures"
OUT.mkdir(exist_ok=True)
df = pd.read_csv(ROOT / "LDSC_master_table.csv")
df = df.rename(columns={"Coefficient_z-score": "Coefficient_z"})

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

mm = df[df.Category == "MM"]
ct_order = (mm.groupby("Celltype")["Coefficient_z"]
              .mean().sort_values(ascending=False).index.tolist())

def z_stars(z):
    if pd.isna(z): return ""
    a = abs(z)
    if a > 3.29: return "***"
    if a > 2.58: return "**"
    if a > 1.96: return "*"
    return ""

vmax_abs = float(np.nanmax(np.abs(df["Coefficient_z"])))
vmin, vmax = -vmax_abs, vmax_abs

fig, axes = plt.subplots(1, 4, figsize=(13.5, 5.2),
                         gridspec_kw={"wspace": 0.45})

cmap = plt.cm.RdBu_r
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

for j, pheno in enumerate(PHENO_ORDER):
    ax = axes[j]
    sub = df[df.Phenotype == pheno]
    mat = np.full((len(ct_order), 2), np.nan)
    for i, ct in enumerate(ct_order):
        for k, mthd in enumerate(METHODS):
            row = sub[(sub.Celltype == ct) & (sub.Category == mthd)]
            if len(row):
                mat[i, k]  = row["Coefficient_z"].iloc[0]

    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)
    for i in range(len(ct_order)):
        for k in range(2):
            if np.isnan(mat[i, k]):
                continue
            txt = f"{mat[i,k]:.2f}{z_stars(mat[i,k])}"
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
cb.set_label("Coefficient z-score", fontsize=10)
cb.ax.tick_params(labelsize=8)

fig.suptitle("LDSC Partitioned h$^2$ Coefficient z-score - Standard vs Suggested (v3, snapMerge DS)",
             fontsize=12, y=1.00)
fig.text(0.5, -0.01,
         "*** |z|>3.29   ** |z|>2.58   * |z|>1.96   "
         "(two-sided normal: p<0.001 / 0.01 / 0.05)",
         ha="center", fontsize=8.5, style="italic")

plt.savefig(OUT / "fig_coef_zscore_heatmap_by_phenotype.pdf",
            bbox_inches="tight", dpi=300)
plt.savefig(OUT / "fig_coef_zscore_heatmap_by_phenotype.png",
            bbox_inches="tight", dpi=200)
print("saved:", OUT / "fig_coef_zscore_heatmap_by_phenotype.pdf")

print("\n=== |z| > 1.96 (p<0.05 two-sided) per phenotype ===")
for pheno in PHENO_ORDER:
    sub = df[df.Phenotype == pheno]
    for m in METHODS:
        s = sub[sub.Category == m]
        n_sig = (s["Coefficient_z"].abs() > 1.96).sum()
        n_pos = ((s["Coefficient_z"] > 1.96)).sum()
        n_neg = ((s["Coefficient_z"] < -1.96)).sum()
        print(f"  {PHENO_LABEL[pheno].splitlines()[0]:<14} {METHOD_LABEL[m]:<10}: "
              f"{n_sig}/{len(s)}  (+{n_pos} / -{n_neg})")
