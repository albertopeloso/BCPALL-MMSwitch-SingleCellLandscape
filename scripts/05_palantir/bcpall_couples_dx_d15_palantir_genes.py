#!/usr/bin/env python
"""Script to plot MAGIC-imputed expression for the final list of genes.

"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Iterable, Optional

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scanpy as sc
import glob
import re

# --- Configuration defaults ---
SAMPLES = {
    "PT06": {"h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt06.h5ad"},
    "PT07": {"h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt07.h5ad"},
    "PT13": {"h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt09.h5ad"},
    "PT14": {"h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt10.h5ad"},
    "PT15": {"h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt11.h5ad"},
}

OUT_H5AD_DIR = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/palantir/processed_palantir_h5ad"
DEFAULT_EMBEDDING_BASIS = "X_draw_graph_fa"
DEFAULT_MAGIC_LAYER = "MAGIC_imputed_data"

DEFAULT_GENES = [
    "MME",      # CD10
    "CD14",     # CD14
    "CD19",     # CD19
    "MS4A1",    # CD20
    "CD22",     # CD22
    "CD24",     # CD24
    "CD33",     # CD33
    "CD34",     # CD34
    "PTPRC",    # CD45
    "CD58",     # CD58
    "CD371",    # CD371 / CLEC12A
]

# CD label shown as subtitle under each gene name in the grid
CD_LABELS = {
    "MME":     "CD10",
    "CD14":    "CD14",
    "CD19":    "CD19",
    "MS4A1":   "CD20",
    "CD22":    "CD22",
    "CD24":    "CD24",
    "CD33":    "CD33",
    "CD34":    "CD34",
    "PTPRC":   "CD45/CD45RA",
    "CD58":    "CD58",
    "CD371":   "CD371",
    "CLEC12A": "CD371",
}

FIGURES_BASE_DIR = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/palantir/figures"
DEFAULT_NEW_GENES_BASE_DIR = os.path.join(FIGURES_BASE_DIR, "new_genes")

LOG = logging.getLogger("newgenes")

# Colormap identified from existing output images
sc.settings.colormap = "Spectral_r"


def _choose_h5ad(sample_cfg: dict, prefer_processed: bool = True) -> Optional[str]:
    """Return path to h5ad to use for this sample.

    If prefer_processed is True and a processed h5ad exists in OUT_H5AD_DIR,
    prefer that file; otherwise fall back to the original path in sample_cfg.
    """
    orig = sample_cfg.get("h5ad")
    if not orig:
        return None
    name = os.path.splitext(os.path.basename(orig))[0]
    proc = os.path.join(OUT_H5AD_DIR, f"{name}_palantir_processed.h5ad")
    if prefer_processed and os.path.exists(proc):
        return proc
    if os.path.exists(orig):
        return orig
    return None


def _has_embedding(ad: anndata.AnnData, basis: str) -> bool:
    return basis in getattr(ad, "obsm", {}) or basis in getattr(ad, "uns", {})


def _resolve_embedding(
    ad: anndata.AnnData,
    embedding: str,
    sample_id: str,
    compute_embedding: bool = False,
) -> Optional[str]:
    """Return a valid embedding key for this AnnData, with fallbacks.
    Returns None if no embedding can be found or computed.
    """
    if _has_embedding(ad, embedding):
        return embedding

    LOG.warning(
        "Requested embedding '%s' not present for %s; trying fallbacks",
        embedding,
        sample_id,
    )
    common = ["X_draw_graph_fa", "X_draw_graph", "X_umap", "X_tsne", "X_pca"]
    for cand in common:
        if _has_embedding(ad, cand):
            LOG.info("Using fallback embedding '%s' for sample %s", cand, sample_id)
            return cand

    # Try any 2D obsm key
    for k, v in getattr(ad, "obsm", {}).items():
        try:
            arr = np.asarray(v)
            if arr.ndim == 2 and arr.shape[1] >= 2:
                LOG.info(
                    "Using fallback embedding '%s' for sample %s", k, sample_id
                )
                return k
        except Exception:
            continue

    # Optionally compute one
    if compute_embedding:
        LOG.info(
            "No embedding found; computing a lightweight embedding for %s", sample_id
        )
        try:
            if "X_pca" not in getattr(ad, "obsm", {}):
                sc.pp.pca(ad)
            sc.pp.neighbors(ad, use_rep="X_pca")
            sc.tl.draw_graph(ad)
            return "X_draw_graph_fa"
        except Exception:
            LOG.exception(
                "Failed to compute fallback embedding for %s", sample_id
            )
            return None

    LOG.warning(
        "No suitable embedding found for %s and --compute-embedding not set; skipping",
        sample_id,
    )
    return None


def _plot_grid(
    ad: anndata.AnnData,
    genes: list,
    outpath: str,
    embedding: str,
    magic_layer: str,
    sample_name: str,
    color_map: str = "Spectral_r",
    ncols: int = 5,
    dpi: int = 300,
) -> None:
    """Plot all genes in a 2-row x 5-col grid in a single image."""

    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"]  = "normal"

    n     = len(genes)
    nrows = int(np.ceil(n / ncols))   # 11 genes / 5 cols = 3 rows; 10 genes = 2 rows
    fig_w = ncols * 4                 # 4 inches per column → 20 inches wide
    fig_h = nrows * 4                 # 4 inches per row   → 8 inches tall

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
    axes = np.array(axes).flatten()

    for i, gene in enumerate(genes):
        ax = axes[i]

        gene_ok  = gene in ad.var_names
        layer_ok = magic_layer in ad.layers

        if not gene_ok:
            ax.set_visible(False)
            LOG.warning(
                "Gene %s not found in var_names for %s; skipping panel",
                gene, sample_name,
            )
            continue

        if not layer_ok:
            ax.set_visible(False)
            LOG.warning(
                "Layer %s not found for %s; skipping panel",
                magic_layer, sample_name,
            )
            continue

        try:
            sc.pl.embedding(
                ad,
                basis=embedding,
                layer=magic_layer,
                color=[gene],
                title="",
                frameon=True,
                wspace=0.4,
                legend_fontweight="normal",
                legend_fontoutline=20,
                color_map=color_map,
                ax=ax,
                show=False,
            )
        except Exception:
            LOG.exception(
                "sc.pl.embedding failed for gene %s sample %s", gene, sample_name
            )
            ax.set_visible(False)
            continue

        # Italic gene name + CD label underneath
        cd_label  = CD_LABELS.get(gene, "")
        title_str = f"$\\it{{{gene}}}$"
        if cd_label:
            title_str += f"\n{cd_label}"
        ax.set_title(title_str, fontsize=14, fontfamily="Karla", pad=6)
        ax.set_xlabel("FDL1", fontsize=10)
        ax.set_ylabel("FDL2", fontsize=10)

    # Hide unused axes (last cell when n=11 and ncols=5 → 3rd row has 1 empty)
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(sample_name, fontsize=18, fontfamily="Karla", y=1.01)
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_single(
    ad: anndata.AnnData,
    gene: str,
    outpath: str,
    embedding: str,
    magic_layer: str,
    color_map: str = "Spectral_r",
    dpi: int = 300,
) -> None:
    """Plot a single gene — used when --no-grid is passed."""
    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"]  = "normal"

    plt.figure(figsize=(6, 5))
    try:
        sc.pl.embedding(
            ad,
            basis=embedding,
            layer=magic_layer,
            color=[gene],
            title="",
            frameon=True,
            wspace=4,
            legend_fontweight="normal",
            legend_fontoutline=20,
            color_map=color_map,
            show=False,
        )
        ax = plt.gca()
        cd_label  = CD_LABELS.get(gene, "")
        title_str = f"$\\it{{{gene}}}$"
        if cd_label:
            title_str += f"\n{cd_label}"
        ax.set_title(title_str, fontfamily="Karla", fontsize=24)
        ax.set_xlabel("FDL1")
        ax.set_ylabel("FDL2")
        plt.tight_layout()
        plt.savefig(outpath, dpi=dpi)
    finally:
        plt.close()


def run(
    genes: Iterable[str],
    samples: Optional[Iterable[str]] = None,
    out_base: str = DEFAULT_NEW_GENES_BASE_DIR,
    prefer_processed: bool = True,
    embedding: str = DEFAULT_EMBEDDING_BASIS,
    magic_layer: str = DEFAULT_MAGIC_LAYER,
    color_map: str = "Spectral_r",
    ncols: int = 5,
    dpi: int = 300,
    compute_embedding: bool = False,
    processed_h5ads: Optional[dict] = None,
    grid_mode: bool = True,
) -> None:
    genes = list(genes)
    samples_list = list(samples) if samples is not None else list(SAMPLES.keys())
    os.makedirs(out_base, exist_ok=True)

    # Build the list of (sample_name, h5ad_path) to process
    if processed_h5ads is not None:
        items = list(processed_h5ads.items())
    else:
        items = []
        for s in samples_list:
            cfg = SAMPLES.get(s)
            if cfg is None:
                LOG.warning("Sample %s not configured; skipping", s)
                continue
            h5 = _choose_h5ad(cfg, prefer_processed=prefer_processed)
            if not h5:
                LOG.warning(
                    "No h5ad found for %s (tried processed and original); skipping", s
                )
                continue
            items.append((s, h5))

    for sample_name, h5 in items:
        LOG.info("Loading %s", h5)
        try:
            ad = anndata.read_h5ad(h5)
            LOG.info(
                "Loaded: %s | Layers: %s | Has MAGIC: %s",
                h5,
                list(ad.layers.keys()),
                magic_layer in ad.layers,
            )
        except Exception as e:
            LOG.exception("Failed to read %s: %s", h5, e)
            continue

        emb = _resolve_embedding(ad, embedding, sample_name, compute_embedding)
        if emb is None:
            continue

        outdir = os.path.join(out_base, sample_name)
        os.makedirs(outdir, exist_ok=True)

        if grid_mode:
            # Single image with all genes in a 2-row x 5-col grid
            outpath = os.path.join(outdir, f"{sample_name}_CD_panel_grid.png")
            try:
                _plot_grid(
                    ad,
                    genes=genes,
                    outpath=outpath,
                    embedding=emb,
                    magic_layer=magic_layer,
                    sample_name=sample_name,
                    color_map=color_map,
                    ncols=ncols,
                    dpi=dpi,
                )
                LOG.info("Saved grid %s", outpath)
            except Exception:
                LOG.exception("Failed to plot grid for %s", sample_name)
        else:
            # Individual image per gene
            for g in genes:
                safe_name = g.replace(os.sep, "_")
                outpath = os.path.join(outdir, f"{sample_name}_{safe_name}.png")
                try:
                    _plot_single(
                        ad,
                        g,
                        outpath,
                        embedding=emb,
                        magic_layer=magic_layer,
                        color_map=color_map,
                        dpi=dpi,
                    )
                    LOG.info("Saved %s", outpath)
                except Exception:
                    LOG.exception("Failed to plot %s for %s", g, sample_name)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Save per-gene MAGIC-imputed embedding plots for samples."
    )
    p.add_argument(
        "--genes",
        nargs="+",
        default=DEFAULT_GENES,
        help="List of genes to plot (default: CD antigen panel)",
    )
    p.add_argument(
        "--samples",
        nargs="*",
        help="Subset of sample keys to process (default: all)",
    )
    p.add_argument(
        "--out-dir",
        default=DEFAULT_NEW_GENES_BASE_DIR,
        help="Base output directory for per-sample folders",
    )
    p.add_argument(
        "--no-prefer-processed",
        dest="prefer_processed",
        action="store_false",
        help="Do not prefer processed h5ad even if present",
    )
    p.add_argument(
        "--embedding",
        default=DEFAULT_EMBEDDING_BASIS,
        help="Embedding basis name to plot (default: X_draw_graph_fa)",
    )
    p.add_argument(
        "--magic-layer",
        default=DEFAULT_MAGIC_LAYER,
        help="Layer name containing MAGIC-imputed data",
    )
    p.add_argument(
        "--color-map",
        default="Spectral_r",
        help="Matplotlib colormap for gene expression (default: Spectral_r)",
    )
    p.add_argument(
        "--ncols",
        type=int,
        default=5,
        help="Number of columns in the grid (default: 5 → 2 rows x 5 cols for 10 genes)",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image DPI",
    )
    p.add_argument(
        "--no-grid",
        dest="grid_mode",
        action="store_false",
        help="Save individual images per gene instead of a grid",
    )
    p.add_argument(
        "--gather-processed",
        action="store_true",
        help="Scan processed h5ad directory and use those files",
    )
    p.add_argument(
        "--processed-dir",
        default=OUT_H5AD_DIR,
        help="Directory containing processed palantir h5ads",
    )
    p.add_argument(
        "--compute-embedding",
        action="store_true",
        help="Compute a lightweight embedding if none are present (may be slow)",
    )
    return p


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    parser = _build_parser()
    args = parser.parse_args()

    processed_h5ads = None
    # If no explicit samples are provided, auto-scan the processed h5ad directory
    # and use those files. Also allow explicit --gather-processed.
    if args.samples is None or args.gather_processed:
        processed_h5ads = {}
        pattern = os.path.join(args.processed_dir, "*_palantir_processed.h5ad")
        for p in glob.glob(pattern):
            bn = os.path.basename(p)
            m  = re.match(r"(.+)_palantir_processed\.h5ad$", bn)
            sample_name = m.group(1) if m else bn
            processed_h5ads[sample_name] = p
        if processed_h5ads:
            LOG.info(
                "Found %d processed h5ads in %s; using these files",
                len(processed_h5ads),
                args.processed_dir,
            )
        else:
            LOG.info("No processed h5ads found in %s", args.processed_dir)

    run(
        genes=args.genes,
        samples=args.samples if args.samples else None,
        out_base=args.out_dir,
        prefer_processed=args.prefer_processed,
        embedding=args.embedding,
        magic_layer=args.magic_layer,
        color_map=args.color_map,
        ncols=args.ncols,
        dpi=args.dpi,
        compute_embedding=args.compute_embedding,
        processed_h5ads=processed_h5ads,
        grid_mode=args.grid_mode,
    )
