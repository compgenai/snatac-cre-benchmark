#!/usr/bin/env python3
"""v3 phase8 .results -> LDSC_master_table.csv (80 rows = 4 trait x 10 cell x 2 method)

MM rows come from v2 symlinks (identical inputs); DS rows are new (standard_snapMerge).
"""
import os
from pathlib import Path
import pandas as pd

# Override:  export LDSC_WORK=/path/to/LDSC/v3
_LDSC_WORK = Path(os.environ.get("LDSC_WORK", "./LDSC_work/v3"))
ROOT  = _LDSC_WORK / "phase8_h2_overlap"
OUT   = _LDSC_WORK / "LDSC_master_table.csv"
CELLS = ["Ventricular_Cardiomyocytes","Atrial_Cardiomyocytes","Fibroblasts",
         "Endothelial","Smooth_Muscle","Macrophages","Nervous_Cells",
         "Myofibroblasts","Primitive_Endoderm","Trophectoderm"]
TRAITS = ["Aragam_2022_CAD_primary","Nauffal_2022_QT",
          "Ntalla_2020_PR_EUR","Roselli_2025_AF_common"]
METHODS = ["DS","MM"]

rows = []
for method in METHODS:
    for trait in TRAITS:
        f = ROOT / f"{method}_{trait}.results"
        df = pd.read_csv(f, sep="\t")
        sel = df[df.Category.str.fullmatch(r"L2_([1-9]|10)")].copy()
        sel = sel.sort_values("Category",
                              key=lambda s: s.str.replace("L2_","").astype(int))
        assert len(sel) == 10, f"{f}: got {len(sel)} L2_ rows"
        sel["Celltype"]  = CELLS
        sel["Category"]  = method
        sel["Phenotype"] = trait
        rows.append(sel)

out = pd.concat(rows, ignore_index=True)
out = out[["Category","Celltype","Phenotype",
           "Prop._SNPs","Prop._h2","Prop._h2_std_error",
           "Enrichment","Enrichment_std_error","Enrichment_p",
           "Coefficient","Coefficient_std_error","Coefficient_z-score"]]
out.to_csv(OUT, index=False)
print(f"wrote {OUT}: {len(out)} rows")
print(out.head())
