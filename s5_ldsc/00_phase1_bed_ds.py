"""Generate phase1 DS BEDs for partitioned_h2_v3 from standard_snapMerge narrowPeak.

Logic (mirrors v2 MM but without AUC/score filter — DS uses ALL peaks):
  1. Read Std_snapMerge_{cell}.narrowPeak
  2. Keep standard chromosomes (chr1-22, chrX, chrY)
  3. Extend +/- 1 kb (start - 1000, end + 1000)
  4. Sort by chrom, start -> write 3-col BED
"""
import os
import pandas as pd
from pathlib import Path

# Overrides:
#   OUT_ROOT   — directory containing standard_snapMerge/*.narrowPeak (S2 · 04 output)
#   LDSC_WORK  — output directory for phase1 hg38 BEDs
SRC_DIR = Path(os.environ.get("OUT_ROOT",  "./output_seung")) / "standard_snapMerge"
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
    src = SRC_DIR / f"Std_snapMerge_{ct}.narrowPeak"
    np_df = pd.read_csv(src, sep="\t", header=None)
    n_orig = len(np_df)

    flt = np_df[np_df[0].isin(STD_CHRS)].copy()
    flt["s"] = (flt[1] - 1000).clip(lower=0)
    flt["e"] = flt[2] + 1000

    bed = flt[[0, "s", "e"]].sort_values([0, "s"]).reset_index(drop=True)
    bed.columns = ["chr", "start", "end"]

    out = OUT_DIR / f"DS_{ct}.bed"
    bed.to_csv(out, sep="\t", header=False, index=False)
    summary_rows.append((ct, n_orig, len(bed)))
    print(f"  {ct:<28} orig={n_orig:>8,}  saved={len(bed):>8,}")

print("\n=== Summary ===")
print(f"{'Celltype':<28} {'orig':>10} {'saved':>10}")
for ct, no, ns in summary_rows:
    print(f"{ct:<28} {no:>10,} {ns:>10,}")
