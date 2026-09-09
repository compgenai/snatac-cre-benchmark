"""Generate phase1 MM BEDs for partitioned_h2_v2 using the corrected
gkmQC AUC>=0.7 + col 8 tie-aware filter + ±1 kb extension.

Logic:
  1. eval.out  → find last topN with AUC > 0.7
  2. top{last_topN}.bed col 5 (1-indexed) min  →  pValue cutoff
  3. Improved_*.narrowPeak col 8 (1-indexed) >= cutoff  →  tie-aware subset
  4. Keep standard chromosomes (chr1-22, chrX, chrY)
  5. Extend ±1 kb (start − 1000, end + 1000)
  6. Sort by chrom, start → write 3-col BED
"""
import os
import pandas as pd
from pathlib import Path

# Overrides:
#   OUT_ROOT   — directory containing gkmQC/MM/ (raw MM centroid narrowPeak)
#   LDSC_WORK  — output directory for phase1 hg38 BEDs consumed by S5
MM_DIR  = Path(os.environ.get("OUT_ROOT",  "./output_seung")) / "gkmQC" / "MM"
OUT_DIR = Path(os.environ.get("LDSC_WORK", "./LDSC_work/v3")) / "phase1_bed_hg38"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CELLS = [
    "Atrial_Cardiomyocytes", "Endothelial", "Fibroblasts", "Macrophages",
    "Myofibroblasts", "Nervous_Cells", "Primitive_Endoderm",
    "Smooth_Muscle", "Trophectoderm", "Ventricular_Cardiomyocytes",
]
STD_CHRS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

summary_rows = []
for ct in CELLS:
    prefix   = f"Improved_{ct}"
    work_dir = MM_DIR / f"{prefix}.gkmqc"

    ed = pd.read_csv(work_dir / f"{prefix}.gkmqc.eval.out", sep="\t",
                     header=None, names=["fa1","fa2","n_seqs","auc","std_err"])
    ed["topN"] = ed["fa1"].str.extract(r"\.top(\d+)\.fa").astype(int)
    last_topN  = int(ed.loc[ed.auc > 0.7, "topN"].max())

    tb     = pd.read_csv(work_dir / f"{prefix}.e300.qc.top{last_topN}.bed",
                         sep="\t", header=None)
    cutoff = float(tb[4].min())

    np_df  = pd.read_csv(MM_DIR / f"{prefix}.narrowPeak", sep="\t", header=None)
    flt    = np_df[(np_df[7] >= cutoff) & (np_df[0].isin(STD_CHRS))].copy()

    flt["s"] = (flt[1] - 1000).clip(lower=0)
    flt["e"] = flt[2] + 1000
    bed = flt[[0, "s", "e"]].sort_values([0, "s"]).reset_index(drop=True)
    bed.columns = ["chr", "start", "end"]

    out = OUT_DIR / f"MM_{ct}.bed"
    bed.to_csv(out, sep="\t", header=False, index=False)
    summary_rows.append((ct, last_topN, cutoff, len(np_df), len(bed)))
    print(f"  {ct:<28} last_topN={last_topN:>3}  "
          f"cutoff={cutoff:>8.4f}  saved={len(bed):>8,}")

print("\n=== Summary ===")
print(f"{'Celltype':<28} {'last_topN':>10} {'cutoff':>10} "
      f"{'orig':>10} {'new_peaks':>10}")
for ct, ln, co, no, npk in summary_rows:
    print(f"{ct:<28} {ln:>10} {co:>10.4f} {no:>10,} {npk:>10,}")
