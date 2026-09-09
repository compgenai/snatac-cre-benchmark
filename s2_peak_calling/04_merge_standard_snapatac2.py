#!/usr/bin/env python
"""Build snapMerge Standard via snapatac2.tl.merge_peaks (half_width=250 -> 501bp).

Replaces the Standard pipeline's bedtools merge step with iterative overlap removal
as implemented in SnapATAC2 v2.8.0. Inputs (per-sample MACS2 narrowPeak) and the
downstream gkmQC pipeline are unchanged.

Per cell type:
  1. Read all per-sample narrowPeak files.
  2. Build polars dict {sample_id: DataFrame} with the schema required by
     snapatac2._snapatac2.py_merge_peaks (chrom/start/end/name/score(u16)/strand/
     signal_value/p_value/q_value/peak).
  3. snap.tl.merge_peaks(peaks, hg38.chrom_sizes, half_width=HALF_WIDTH) ->
     polars.DataFrame with "Peaks" string column ("chr:s-e", 501bp fixed-width).
  4. For each merged interval, map back to the strongest contributing source peak
     (max p_value among source peaks whose summit lies in the interval) and emit
     a 10-col narrowPeak row with that source peak's score/signal/p/q and
     peak_offset = summit - merged_start. This keeps the score metric identical
     to Suggested (-log10 pValue), enabling apples-to-apples comparison.
  5. Write Std_snapMerge_{ct}.narrowPeak (sorted chrom, start).

Trophectoderm has only 8 samples (sample8 is missing on disk); processed as-is.
"""
from __future__ import annotations
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl


# Overrides:
#   DATA_ROOT — directory containing peakcalling_default/ (S2 · 02 output)
#   OUT_ROOT  — directory to write standard_snapMerge/*.narrowPeak
PEAKCALL_ROOT = Path(os.environ.get("DATA_ROOT", "./output_so")) / "peakcalling_default"
OUT_DIR       = Path(os.environ.get("OUT_ROOT",  "./output_seung")) / "standard_snapMerge"
OUT_DIR.mkdir(parents=True, exist_ok=True)
HALF_WIDTH = 250  # 501 bp fixed width (Corces 2018 / SnapATAC2 / ArchR canonical)

CELLTYPES = [
    "Atrial_Cardiomyocytes", "Endothelial", "Fibroblasts", "Macrophages",
    "Myofibroblasts", "Nervous_Cells", "Primitive_Endoderm", "Smooth_Muscle",
    "Trophectoderm", "Ventricular_Cardiomyocytes",
]


def sample_narrowpeak_path(ct: str, s: int) -> Path:
    """Per-sample MACS2 narrowPeak path; returns Path even if missing (caller checks)."""
    return (
        PEAKCALL_ROOT / f"mcluster{ct}.sample{s}" /
        f"macs2_mcluster{ct}.sample{s}" /
        f"mcluster{ct}.sample{s}_peaks.narrowPeak"
    )


def list_sample_files(ct: str) -> list[tuple[str, Path]]:
    """Return [(sample_id, path)] for samples whose narrowPeak exists AND is non-empty.

    Trophectoderm in particular has multiple empty narrowPeak files (sample5/6/9 = 0 lines)
    and sample8 is missing entirely; these are skipped with a warning.
    """
    out = []
    for s in range(1, 10):
        p = sample_narrowpeak_path(ct, s)
        if not p.exists():
            print(f"[{ct}] WARN sample{s}: file missing", flush=True)
            continue
        if p.stat().st_size == 0:
            print(f"[{ct}] WARN sample{s}: empty file (skipped)", flush=True)
            continue
        out.append((f"sample{s}", p))
    return out


def load_narrowpeak_polars(path: Path) -> pl.DataFrame:
    """Read MACS2 narrowPeak (10 cols) into a polars DataFrame.

    MACS2 occasionally emits col5 (score) > 1000 (e.g. 133474), so it is read
    as UInt32 here and clipped to UInt16 only when building the merge_peaks input
    (whose Rust schema requires UInt16). The original wide score is preserved
    in the output narrowPeak.
    """
    df = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom", "start", "end", "name", "score", "strand",
            "signal_value", "p_value", "q_value", "peak",
        ],
        schema_overrides={
            "chrom": pl.Utf8, "start": pl.UInt64, "end": pl.UInt64,
            "name": pl.Utf8, "score": pl.UInt32, "strand": pl.Utf8,
            "signal_value": pl.Float64, "p_value": pl.Float64,
            "q_value": pl.Float64, "peak": pl.UInt64,
        },
    )
    return df


def build_celltype(ct: str, chrom_sizes: dict[str, int]) -> Path:
    """Run merge_peaks for one cell type and write Std_snapMerge_{ct}.narrowPeak.

    Returns the output path.
    """
    import snapatac2 as snap

    sample_files = list_sample_files(ct)
    if not sample_files:
        raise RuntimeError(f"[{ct}] no per-sample narrowPeak found")

    print(f"[{ct}] {len(sample_files)} samples:", [s for s, _ in sample_files], flush=True)

    t0 = time.time()
    per_sample: dict[str, pl.DataFrame] = {}
    src_rows: list[pl.DataFrame] = []
    # snapatac2 v2.8.0 advertises UInt16 for score but internally casts to u8 (0-255).
    # Values >255 trigger a Rust TryFromInt panic; clip here. Score is ignored by the
    # merge_peaks algorithm (ranking is by p_value), and the original wide score is
    # preserved separately in src_rows for the post-merge mapping step.
    U8_MAX = 255
    for sid, p in sample_files:
        df = load_narrowpeak_polars(p)
        per_sample[sid] = df.with_columns(
            pl.col("score").clip(upper_bound=U8_MAX).cast(pl.UInt16)
        )
        # Keep a tagged copy (with original wide score) for the source mapping step.
        src_rows.append(
            df.with_columns(
                pl.lit(sid).alias("sample"),
                (pl.col("start") + pl.col("peak")).alias("summit"),
            )
        )
    src_all = pl.concat(src_rows)
    t_read = time.time() - t0
    print(f"[{ct}] read {len(src_all):,} source peaks in {t_read:.1f}s", flush=True)

    # --- merge ---
    t1 = time.time()
    merged = snap.tl.merge_peaks(per_sample, chrom_sizes, half_width=HALF_WIDTH)
    t_merge = time.time() - t1
    print(f"[{ct}] merged -> {len(merged):,} intervals in {t_merge:.1f}s", flush=True)

    # --- parse "chr:s-e" -> chrom/start/end ---
    parsed = merged.with_columns(
        pl.col("Peaks").str.extract(r"^(.+):", 1).alias("chrom"),
        pl.col("Peaks").str.extract(r":(\d+)-", 1).cast(pl.Int64).alias("m_start"),
        pl.col("Peaks").str.extract(r"-(\d+)$", 1).cast(pl.Int64).alias("m_end"),
    ).select(["chrom", "m_start", "m_end"])

    # --- map back: for each merged interval find max-p_value source peak whose summit
    # lies within [m_start, m_end). Use polars join_asof + filter for speed.
    # Strategy: for each source row, find the merged interval containing its summit
    # (per-chrom sorted), then group_by interval to pick max p_value.
    t2 = time.time()

    merged_sorted = parsed.sort(["chrom", "m_start"]).with_row_index("mid")
    # Per-chrom asof join: join each source row with the latest m_start <= summit.
    src_sorted = src_all.select(["chrom", "summit", "score", "signal_value",
                                 "p_value", "q_value", "sample"]).sort(["chrom", "summit"])

    joined = src_sorted.join_asof(
        merged_sorted.sort(["chrom", "m_start"]),
        left_on="summit", right_on="m_start", by="chrom", strategy="backward",
    )
    # Keep only sources whose summit falls strictly within the matched interval.
    joined = joined.filter(
        pl.col("mid").is_not_null() &
        (pl.col("summit") >= pl.col("m_start")) &
        (pl.col("summit") < pl.col("m_end"))
    )

    # Pick the source with max p_value per merged interval (tie-break: lower summit).
    best = (
        joined.sort(["mid", "p_value", "summit"], descending=[False, True, False])
        .group_by("mid", maintain_order=True)
        .first()
    )

    # Build 10-col narrowPeak rows.
    # peak_offset = summit - m_start (clipped to >=0 for edge-clipped intervals).
    out_df = best.with_columns(
        (pl.col("summit") - pl.col("m_start")).clip(lower_bound=0).alias("peak_offset"),
    ).select([
        pl.col("chrom"),
        pl.col("m_start").cast(pl.Int64).alias("start"),
        pl.col("m_end").cast(pl.Int64).alias("end"),
        (pl.lit(f"snapMerge_{ct}_") + pl.col("mid").cast(pl.Utf8)).alias("name"),
        pl.col("score").cast(pl.UInt32).alias("score"),
        pl.lit(".").alias("strand"),
        pl.col("signal_value"),
        pl.col("p_value"),
        pl.col("q_value"),
        pl.col("peak_offset").cast(pl.Int64).alias("peak"),
    ]).sort(["chrom", "start"])

    n_lost = len(merged) - len(out_df)
    if n_lost > 0:
        # Some merged intervals had no source summit landing inside (should be 0
        # under merge_peaks semantics; sanity check only).
        print(f"[{ct}] WARN: {n_lost} merged intervals lacked source-summit mapping", flush=True)

    out_path = OUT_DIR / f"Std_snapMerge_{ct}.narrowPeak"
    out_df.write_csv(out_path, separator="\t", include_header=False)
    t_map = time.time() - t2
    print(f"[{ct}] mapped + wrote {len(out_df):,} rows in {t_map:.1f}s -> {out_path}", flush=True)
    print(f"[{ct}] total wall {time.time() - t0:.1f}s", flush=True)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--celltype", "-c", default="all",
        help="cell type to process, or 'all' (default)",
    )
    args = ap.parse_args()

    import snapatac2 as snap
    chrom_sizes = snap.genome.hg38.chrom_sizes

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.celltype == "all":
        cts = CELLTYPES
    else:
        if args.celltype not in CELLTYPES:
            raise SystemExit(f"unknown cell type: {args.celltype}")
        cts = [args.celltype]

    for ct in cts:
        try:
            build_celltype(ct, chrom_sizes)
        except Exception as e:
            print(f"[{ct}] FAILED: {e}", file=sys.stderr, flush=True)
            raise


if __name__ == "__main__":
    main()
