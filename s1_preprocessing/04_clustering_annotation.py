#!/usr/bin/env python3
"""
S1 · 04 — Merge per-sample AnnData, cluster, annotate cell types, export fragments.

Runs the full end-to-end clustering + annotation + fragment-export pipeline:
    1. Merge per-sample h5ad files into a single AnnDataSet.
    2. select_features(n_features=200000) → snap.tl.spectral (30-d).
    3. snap.pp.harmony(batch="sample", use_dims=30, max_iter_harmony=20).
    4. snap.pp.knn(k=10) → snap.tl.leiden(random_state=0) → snap.tl.umap.
    5. Build a gene-matrix, render the marker-gene dotplot, and (optionally)
       write UMAP/dotplot PDFs.
    6. Map Leiden clusters → cell types via ``--celltype-mapping`` JSON.
    7. Write per-sample updated h5ads (with obs["celltype"], obs["leiden"],
       obs["sample"], obs["sample_cluster"]).
    8. Merge them into a new AnnDataSet and export per-``sample_cluster``
       fragment BEDs under ``$DATA_ROOT/fragments/mcluster<CELLTYPE>.<SAMPLE>_fragments.bed``
       (consumed by ``s2_peak_calling/01_make_tn5_bed.sh``).

Environment overrides:
    PROJECT_ROOT   project root (default ".")
    DATA_ROOT      "$PROJECT_ROOT/output_so" by default

The Leiden → cell type mapping is supplied via ``--celltype-mapping``
(default: ``celltype_mapping.json`` next to this script). Marker genes for the
dotplot come from ``--marker-genes`` (default: ``marker_genes.json``).
"""
from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

# Limit thread pools before importing numeric libs (matches the notebook).
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("RAYON_NUM_THREADS", "1")
os.environ.setdefault("MALLOC_ARENA_MAX", "2")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import snapatac2 as snap


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    project_root = Path(os.environ.get("PROJECT_ROOT", "."))
    data_root = Path(os.environ.get("DATA_ROOT", project_root / "output_so"))

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", type=Path, default=data_root)
    p.add_argument("--samples", nargs="+", default=[f"sample{i}" for i in range(1, 10)])
    p.add_argument("--aggr-dir", type=Path, default=None,
                   help="Per-sample h5ad directory (default $DATA_ROOT/aggr)")
    p.add_argument("--aggr-update-dir", type=Path, default=None,
                   help="Output directory for per-sample updated h5ads and merged AnnDataSet "
                        "(default $DATA_ROOT/aggr_update)")
    p.add_argument("--fragments-dir", type=Path, default=None,
                   help="Output directory for per-cluster fragment BEDs (default $DATA_ROOT/fragments)")
    p.add_argument("--figures-dir", type=Path, default=None,
                   help="Output directory for UMAP / dotplot PDFs (default $DATA_ROOT/qc_figures/clustering)")
    p.add_argument("--celltype-mapping", type=Path, default=SCRIPT_DIR / "celltype_mapping.json",
                   help="JSON mapping Leiden cluster id (string) → celltype name")
    p.add_argument("--marker-genes", type=Path, default=SCRIPT_DIR / "marker_genes.json",
                   help="JSON marker-gene groups for the dotplot")
    p.add_argument("--n-features", type=int, default=200_000)
    p.add_argument("--use-dims", type=int, default=30)
    p.add_argument("--knn-k", type=int, default=10)
    p.add_argument("--max-iter-harmony", type=int, default=20)
    p.add_argument("--random-state", type=int, default=0)
    return p.parse_args()


def load_json(path: Path) -> dict:
    with open(path) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}


def merge_h5ads(aggr_dir: Path, samples: list[str], merged_path: Path):
    file_paths = [(s, str(aggr_dir / f"{s}.h5ad")) for s in samples]
    for _, p in file_paths:
        if not Path(p).exists():
            raise FileNotFoundError(f"Missing per-sample h5ad: {p}")
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    logging.info("Merging %d per-sample h5ads → %s", len(file_paths), merged_path)
    return snap.AnnDataSet(adatas=file_paths, filename=str(merged_path))


def cluster(data, args: argparse.Namespace) -> None:
    logging.info("select_features n_features=%d", args.n_features)
    snap.pp.select_features(data, n_features=args.n_features)
    logging.info("spectral embedding")
    snap.tl.spectral(data)
    logging.info("harmony batch='sample' use_dims=%d max_iter=%d", args.use_dims, args.max_iter_harmony)
    snap.pp.harmony(data, batch="sample", use_dims=args.use_dims, max_iter_harmony=args.max_iter_harmony)
    logging.info("kNN k=%d", args.knn_k)
    snap.pp.knn(data, n_neighbors=args.knn_k, use_dims=args.use_dims, use_rep="X_spectral_harmony")
    logging.info("Leiden random_state=%d", args.random_state)
    snap.tl.leiden(data, random_state=args.random_state)
    logging.info("UMAP")
    snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=args.use_dims)


def render_umaps(gene_matrix, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for color, name in (("leiden", "umap_leiden"), ("sample", "umap_sample")):
        sc.pl.umap(gene_matrix, use_raw=False, color=color, show=False)
        out = out_dir / f"{name}.pdf"
        plt.savefig(out, bbox_inches="tight")
        plt.close("all")
        logging.info("Wrote %s", out)


def render_dotplot(gene_matrix, marker_genes: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    dotplot = sc.pl.dotplot(
        gene_matrix,
        var_names=marker_genes,
        groupby="leiden",
        standard_scale="var",
        dot_min=0.1,
        dot_max=1,
        cmap="Blues",
        figsize=(12, 10),
        show=False,
    )
    fig = getattr(dotplot, "figure", None) or dotplot["dotplot_ax"].figure
    out = out_dir / "marker_dotplot.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    logging.info("Wrote %s", out)


def update_persample_h5ads(gene_matrix, aggr_dir: Path, aggr_update_dir: Path,
                            samples: list[str], mapping: dict) -> list[tuple[str, str]]:
    gene_matrix = gene_matrix[~gene_matrix.obs.index.duplicated(keep="first")].copy()
    leiden = gene_matrix.obs["leiden"].astype(str)
    sample = gene_matrix.obs["sample"].astype(str)
    cell_types = leiden.map(lambda c: mapping.get(c, "Unknown"))
    barcodes = gene_matrix.obs.index.astype(str)
    df = pd.DataFrame({"leiden": leiden, "sample": sample, "celltype": cell_types}, index=barcodes)

    aggr_update_dir.mkdir(parents=True, exist_ok=True)
    updated_files: list[tuple[str, str]] = []
    for s in samples:
        src = aggr_dir / f"{s}.h5ad"
        if not src.exists():
            logging.warning("skipping %s (missing %s)", s, src)
            continue
        ad = sc.read_h5ad(src)
        ad = ad[~ad.obs.index.duplicated(keep="first")].copy()
        merged = df.reindex(ad.obs.index.astype(str))
        merged["celltype"] = merged["celltype"].fillna("Unknown").astype(str)
        merged["leiden"] = merged["leiden"].fillna("NA").astype(str)
        merged["sample"] = merged["sample"].fillna(s).astype(str)

        ad.obs["celltype"] = pd.Categorical(merged["celltype"].values)
        ad.obs["leiden"] = pd.Categorical(merged["leiden"].values)
        ad.obs["sample"] = pd.Categorical(merged["sample"].values)
        ad.obs["sample_cluster"] = ad.obs["celltype"].astype(str) + "." + ad.obs["sample"].astype(str)

        out_path = aggr_update_dir / f"{s}_updated.h5ad"
        ad.write(out_path)
        updated_files.append((s, str(out_path)))
        logging.info("Wrote %s (n=%d)", out_path, ad.n_obs)
    return updated_files


def export_fragments(updated_files: list[tuple[str, str]], merged_out: Path, fragments_dir: Path) -> None:
    data = snap.AnnDataSet(adatas=updated_files, filename=str(merged_out))

    labels = []
    for _, p in updated_files:
        ad = sc.read_h5ad(p)
        if "sample_cluster" not in ad.obs:
            raise ValueError(f"{p} missing 'sample_cluster'")
        labels.append(ad.obs["sample_cluster"].astype(str))
    groupby_series = pd.concat(labels)

    fragments_dir.mkdir(parents=True, exist_ok=True)
    logging.info("Exporting per-(sample × celltype) fragments → %s", fragments_dir)
    snap.ex.export_fragments(
        data,
        groupby=groupby_series,
        out_dir=str(fragments_dir),
        prefix="mcluster",
        suffix="_fragments.bed",
    )


def main() -> None:
    logging.basicConfig(format="[%(levelname)s] %(message)s", level=logging.INFO)
    args = parse_args()

    aggr_dir = args.aggr_dir or (args.data_root / "aggr")
    aggr_update_dir = args.aggr_update_dir or (args.data_root / "aggr_update")
    fragments_dir = args.fragments_dir or (args.data_root / "fragments")
    figures_dir = args.figures_dir or (args.data_root / "qc_figures" / "clustering")

    mapping = load_json(args.celltype_mapping)
    marker_genes = load_json(args.marker_genes)

    merged_path = aggr_dir / "colon.h5ads"
    data = merge_h5ads(aggr_dir, args.samples, merged_path)

    cluster(data, args)

    logging.info("Building gene matrix")
    gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)
    gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]
    gene_matrix.obs["leiden"] = data.obs["leiden"]

    render_umaps(gene_matrix, figures_dir)
    render_dotplot(gene_matrix, marker_genes, figures_dir)

    data.close()

    updated_files = update_persample_h5ads(gene_matrix, aggr_dir, aggr_update_dir,
                                            args.samples, mapping)
    if not updated_files:
        logging.error("No updated per-sample files were written; aborting fragment export.")
        return

    export_fragments(updated_files, aggr_update_dir / "colon_updated.h5ads", fragments_dir)
    logging.info("Done.")


if __name__ == "__main__":
    main()
