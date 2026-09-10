#!/usr/bin/env python3
"""
S1 · 02 — Per-sample QC and cell filtering.

Reads Cell Ranger ATAC per-sample fragment files, computes TSS enrichment,
and writes a filtered barcode table per sample. Also caches a per-sample
h5ad (raw import + TSSE) for the doublet-removal step.

Filter (applied per sample, uniform thresholds):
    is__cell_barcode == 1
    TSS_enrichment >= TSSE_MIN                (default 6)
    passed_filters + 1 in [UMI_MIN, UMI_MAX]  (default 5,000 – 50,000)

Optional QC figures (TSSE violin, fragment-size distribution, mitochondrial
fraction) are written to ``<qc_fig_dir>`` when ``--qc-figures`` is passed.

Environment overrides:
    PROJECT_ROOT   project root (default ".")
    DATA_ROOT      "$PROJECT_ROOT/output_so" by default
"""
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import snapatac2 as snap


DEFAULT_TSSE_MIN = 6.0
DEFAULT_UMI_MIN = 5_000
DEFAULT_UMI_MAX = 50_000


def parse_args() -> argparse.Namespace:
    project_root = Path(os.environ.get("PROJECT_ROOT", "."))
    data_root = Path(os.environ.get("DATA_ROOT", project_root / "output_so"))

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", type=Path, default=data_root,
                   help="Root directory holding cell_ranger_output/ (default $DATA_ROOT)")
    p.add_argument("--samples", nargs="+", default=[f"sample{i}" for i in range(1, 10)],
                   help="Sample IDs (folder names under cell_ranger_output/)")
    p.add_argument("--tsse-min", type=float, default=DEFAULT_TSSE_MIN)
    p.add_argument("--umi-min", type=int, default=DEFAULT_UMI_MIN)
    p.add_argument("--umi-max", type=int, default=DEFAULT_UMI_MAX)
    p.add_argument("--h5ad-dir", type=Path, default=None,
                   help="Output directory for per-sample cached h5ad (default: $DATA_ROOT/h5ad)")
    p.add_argument("--filtered-dir", type=Path, default=None,
                   help="Output directory for per-sample filtered CSV (default: $DATA_ROOT/filtered_samples)")
    p.add_argument("--qc-figures", action="store_true",
                   help="Also render per-sample QC figures (TSSE violin, fragment-size, mito)")
    p.add_argument("--qc-fig-dir", type=Path, default=None,
                   help="Directory for QC figures (default: $DATA_ROOT/qc_figures)")
    return p.parse_args()


def fragment_file(data_root: Path, sample: str) -> Path:
    return data_root / "cell_ranger_output" / sample / "outs" / "fragments.tsv.gz"


def singlecell_csv(data_root: Path, sample: str) -> Path:
    return data_root / "cell_ranger_output" / sample / "outs" / "singlecell.csv"


def import_with_tsse(frag_path: Path):
    """Import a fragment file and compute TSS enrichment (no on-disk h5ad backing)."""
    adata = snap.pp.import_fragments(
        str(frag_path),
        chrom_sizes=snap.genome.hg38,
        sorted_by_barcode=False,
    )
    snap.metrics.tsse(adata, snap.genome.hg38)
    return adata


def cache_h5ad(frag_path: Path, out_h5ad: Path):
    """Import a fragment file, compute TSSE, and persist to h5ad."""
    if out_h5ad.exists():
        out_h5ad.unlink()
    adata = import_with_tsse(frag_path)
    adata.write_h5ad(str(out_h5ad))
    return adata


def render_fragment_size_grid(adatas: dict, out_png: Path):
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    for i, (sample, adata) in enumerate(adatas.items()):
        ax = axes[i // 3, i % 3]
        snap.metrics.frag_size_distr(adata, add_key="frag_size_distr", max_recorded_size=1000)
        dist = adata.uns["frag_size_distr"]
        ax.plot(range(1, len(dist)), dist[1:])
        ax.set_xlabel("Fragment size")
        ax.set_ylabel("Count")
        ax.set_title(sample)
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_tsse_violin(tsse_scores: dict, tsse_min: float, out_png: Path):
    labels = sorted(tsse_scores.keys())
    scores = [tsse_scores[k] for k in labels]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.violinplot(scores, positions=x, showmeans=False, showmedians=True)
    ax.axhline(y=tsse_min, color="red", linestyle="--", label=f"TSSe = {tsse_min}")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30)
    ax.set_ylabel("TSS Enrichment Score")
    ax.set_title("TSSe distribution across samples")
    ax.legend()
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_tsse_individual(adatas: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = []
    for sample, adata in adatas.items():
        plt.close("all")
        snap.pl.tsse(adata, interactive=False)
        fig = plt.gcf()
        out_png = out_dir / f"{sample}_tsse.png"
        fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        saved.append(out_png)
    return saved


def render_mito_scatter(qc_dfs: dict, out_png: Path, umi_col: str = "passed_filters",
                        mito_col: str = "mitochondrial", y_cutoff: float = 0.1):
    samples = sorted(qc_dfs.keys())
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    for i, sample in enumerate(samples):
        ax = axes[i // 3, i % 3]
        qc = qc_dfs[sample]
        cells = qc.query("is__cell_barcode == 1")
        for subset, label, size, alpha in ((qc, "All", 1, 0.4), (cells, "Filtered", 1, 1.0)):
            frag = subset[mito_col] + 1
            umi = subset[umi_col] + 1
            frac = frag / umi
            logumi = np.log10(umi + 1)
            ax.scatter(logumi, frac, s=size, alpha=alpha, label=label)
        ax.set_xlim(2, 5)
        ax.set_ylim(0, 1)
        ax.set_title(sample)
        ax.set_xlabel("log10(UMI+1)")
        ax.set_ylabel(f"{mito_col}/(UMI + {mito_col})")
        ax.axhline(y=y_cutoff, color="r", linestyle="-", lw=0.5)
        ax.axvline(x=3, color="r", linestyle="-", lw=0.5)
        if i == 0:
            ax.legend(markerscale=6, fontsize=8)
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)
    args = parse_args()

    h5ad_dir = args.h5ad_dir or (args.data_root / "h5ad")
    filtered_dir = args.filtered_dir or (args.data_root / "filtered_samples")
    qc_fig_dir = args.qc_fig_dir or (args.data_root / "qc_figures")

    h5ad_dir.mkdir(parents=True, exist_ok=True)
    filtered_dir.mkdir(parents=True, exist_ok=True)

    tsse_scores: dict[str, list[float]] = {}
    qc_dfs: dict[str, pd.DataFrame] = {}
    adatas_for_figures: dict = {}

    summary_rows = []
    for sample in args.samples:
        frag = fragment_file(args.data_root, sample)
        qc = singlecell_csv(args.data_root, sample)
        out_h5ad = h5ad_dir / f"{sample}.h5ad"

        logging.info("[%s] importing fragments and computing TSSe", sample)
        adata = cache_h5ad(frag, out_h5ad)
        tsse_scores[sample] = list(adata.obs["tsse"])
        if args.qc_figures:
            adatas_for_figures[sample] = adata

        logging.info("[%s] loading singlecell.csv", sample)
        qc_df = pd.read_csv(qc, index_col="barcode")
        qc_df["TSS_enrichment"] = qc_df.index.map(dict(zip(adata.obs_names, adata.obs["tsse"])))
        qc_dfs[sample] = qc_df

        cells = qc_df.query("is__cell_barcode == 1")
        filtered = cells[
            (cells["TSS_enrichment"] >= args.tsse_min)
            & (cells["passed_filters"] + 1 >= args.umi_min)
            & (cells["passed_filters"] + 1 <= args.umi_max)
        ]
        out_csv = filtered_dir / f"{sample}_filtered.csv"
        filtered.to_csv(out_csv, index=True)
        pre, post = len(qc_df), len(filtered)
        retention = post / pre if pre else 0.0
        summary_rows.append((sample, pre, post, retention))
        logging.info("[%s] filtered %d/%d cells (retention %.4f) → %s", sample, post, pre, retention, out_csv)

    summary = pd.DataFrame(summary_rows, columns=["sample", "pre_cells", "post_cells", "retention"])
    logging.info("\n%s", summary.to_string(index=False))
    logging.info("TOTAL pre=%d post=%d retention=%.4f",
                 summary["pre_cells"].sum(),
                 summary["post_cells"].sum(),
                 summary["post_cells"].sum() / summary["pre_cells"].sum())

    if args.qc_figures:
        qc_fig_dir.mkdir(parents=True, exist_ok=True)
        logging.info("Rendering QC figures → %s", qc_fig_dir)
        render_fragment_size_grid(adatas_for_figures, qc_fig_dir / "fragment_size_grid.png")
        render_tsse_violin(tsse_scores, args.tsse_min, qc_fig_dir / "tsse_violin.png")
        render_tsse_individual(adatas_for_figures, qc_fig_dir / "tsse_per_sample")
        render_mito_scatter(qc_dfs, qc_fig_dir / "mito_scatter.png")


if __name__ == "__main__":
    main()
