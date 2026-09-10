#!/usr/bin/env python3
"""
S1 · 03 — Doublet removal per sample.

For each sample, reads the per-sample filtered CSV from S1 · 02, imports the
Cell Ranger fragments backed by an on-disk h5ad, subsets to the filtered
barcodes, builds a 500-bp tile matrix, selects features (n=250,000), runs
``snap.pp.scrublet`` with default parameters, and filters doublets at
probability > 0.5 (SnapATAC2 default).

Output: one AnnData per sample at ``$DATA_ROOT/cell_qc/doublets/<sample>_doublets.h5ad``.

Environment overrides:
    PROJECT_ROOT   project root (default ".")
    DATA_ROOT      "$PROJECT_ROOT/output_so" by default
"""
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import pandas as pd
import snapatac2 as snap


def parse_args() -> argparse.Namespace:
    project_root = Path(os.environ.get("PROJECT_ROOT", "."))
    data_root = Path(os.environ.get("DATA_ROOT", project_root / "output_so"))

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", type=Path, default=data_root,
                   help="Root directory holding cell_ranger_output/ and filtered_samples/ (default $DATA_ROOT)")
    p.add_argument("--samples", nargs="+", default=[f"sample{i}" for i in range(1, 10)],
                   help="Sample IDs")
    p.add_argument("--filtered-dir", type=Path, default=None,
                   help="Per-sample filtered CSV directory (default $DATA_ROOT/filtered_samples)")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Output directory for doublet-removed h5ads (default $DATA_ROOT/cell_qc/doublets)")
    p.add_argument("--n-features", type=int, default=250_000,
                   help="snap.pp.select_features n_features (default 250,000)")
    return p.parse_args()


def fragment_file(data_root: Path, sample: str) -> Path:
    return data_root / "cell_ranger_output" / sample / "outs" / "fragments.tsv.gz"


def process_sample(sample: str, data_root: Path, filtered_dir: Path, out_dir: Path,
                   n_features: int) -> int | None:
    filtered_csv = filtered_dir / f"{sample}_filtered.csv"
    if not filtered_csv.exists():
        logging.warning("[%s] filtered CSV missing (%s), skipping", sample, filtered_csv)
        return None

    frag = fragment_file(data_root, sample)
    if not frag.exists():
        logging.warning("[%s] fragment file missing (%s), skipping", sample, frag)
        return None

    filtered = pd.read_csv(filtered_csv, index_col="barcode")
    kept_barcodes = set(filtered.index)

    out_h5ad = out_dir / f"{sample}_doublets.h5ad"
    logging.info("[%s] importing fragments → %s", sample, out_h5ad)
    data = snap.pp.import_data(
        str(frag),
        chrom_sizes=snap.genome.hg38,
        file=str(out_h5ad),
        sorted_by_barcode=False,
    )

    mask = [bc in kept_barcodes for bc in data.obs_names]
    data.subset(obs_indices=mask)

    logging.info("[%s] add_tile_matrix (500 bp)", sample)
    snap.pp.add_tile_matrix(data)
    logging.info("[%s] select_features n_features=%d", sample, n_features)
    snap.pp.select_features(data, n_features=n_features)
    logging.info("[%s] scrublet + filter_doublets", sample)
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)

    logging.info("[%s] cells after doublet removal: %d", sample, data.shape[0])
    return data.shape[0]


def main() -> None:
    logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)
    args = parse_args()

    filtered_dir = args.filtered_dir or (args.data_root / "filtered_samples")
    out_dir = args.out_dir or (args.data_root / "cell_qc" / "doublets")
    out_dir.mkdir(parents=True, exist_ok=True)

    totals = {}
    for sample in args.samples:
        n = process_sample(sample, args.data_root, filtered_dir, out_dir, args.n_features)
        if n is not None:
            totals[sample] = n

    if totals:
        logging.info("Summary (cells after doublet removal):")
        for sample, n in sorted(totals.items()):
            logging.info("  %s: %d", sample, n)
        logging.info("  TOTAL: %d", sum(totals.values()))


if __name__ == "__main__":
    main()
