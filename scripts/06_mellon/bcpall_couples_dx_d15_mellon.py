#!/usr/bin/env python
# coding: utf-8

"""
Mellon density analysis on multiple samples
- save-only plotting (no GUI)
- batch processing of multiple h5ad files
- uses Karla font and BCPALL patient names
- uses Palantir-processed h5ad files to reproduce previous log_density ranges
"""

# =========================
# Imports & global setup
# =========================

import scanpy as sc
import numpy as np
import mellon
import palantir
import os
import logging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---- Set global font to Karla ----
rcParams["font.family"] = "Karla"

# ---- Logging ----
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

# =========================
# Configuration
# =========================

SAMPLES = {
    "PT06": {"h5ad": "PT06_palantir_processed.h5ad"},
    "PT07": {"h5ad": "PT07_palantir_processed.h5ad"},
    "PT09": {"h5ad": "PT09_palantir_processed.h5ad"},
    "PT10": {"h5ad": "PT10_palantir_processed.h5ad"},
    "PT11": {"h5ad": "PT11_palantir_processed.h5ad"},
}

PALANTIR_DIR = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/palantir/processed_palantir_h5ad"
OUT_H5AD_DIR = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/mellon/processed_h5ad"
FIGURES_BASE_DIR = "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/mellon/figures"

# Mapping PT → BCPALL codes
pt_to_bcpall = {
    "PT06": "BCPALL#06",
    "PT07": "BCPALL#07",
    "PT09": "BCPALL#09",
    "PT10": "BCPALL#10",
    "PT11": "BCPALL#11",
}

# =========================
# Pipeline
# =========================

def run_mellon_pipeline(sample_name: str, cfg: dict):

    h5ad_path = os.path.join(PALANTIR_DIR, cfg["h5ad"])

    # Output dirs
    figdir = os.path.join(FIGURES_BASE_DIR, sample_name)
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    os.makedirs(OUT_H5AD_DIR, exist_ok=True)

    logger.info(f"=== Running Mellon for {sample_name} ===")

    # ---------------------
    # Load Palantir-processed AnnData
    # ---------------------
    ad = sc.read(h5ad_path)

    # ---------------------
    # Ensure PCA exists
    # ---------------------
    if "X_pca" not in ad.obsm:
        logger.info("PCA missing — computing PCA")
        sc.pp.pca(ad)

    # ---------------------
    # Ensure diffusion maps exist
    # ---------------------
    if "DM_EigenVectors" not in ad.obsm:
        logger.info("Diffusion maps missing — computing them")
        palantir.utils.run_diffusion_maps(ad, pca_key="X_pca", n_components=10)
    embedding = ad.obsm["DM_EigenVectors"]

    # ---------------------
    # Run Mellon density
    # ---------------------
    logger.info("Running Mellon density estimation")
    model = mellon.DensityEstimator()
    log_density = model.fit_predict(embedding)

    ad.obs["mellon_log_density"] = log_density
    ad.obs["mellon_log_density_clipped"] = np.clip(
        log_density, *np.quantile(log_density, [0.05, 1])
    )
    ad.uns["mellon_log_density_function"] = model.predict.to_dict()

    # ---------------------
    # Ensure ForceAtlas layout exists
    # ---------------------
    if "X_draw_graph_fa" not in ad.obsm:
        logger.info("ForceAtlas layout missing — computing graph layout")
        sc.pp.neighbors(ad)
        sc.tl.draw_graph(ad)

    # ---------------------
    # Plot Mellon densities
    # ---------------------
    sc.set_figure_params(dpi_save=300)

    for key, suffix in [("mellon_log_density", "mellon_density"),
                        ("mellon_log_density_clipped", "mellon_density_clipped")]:
        sc.pl.embedding(
            ad, basis="X_draw_graph_fa", color=key, frameon=True, show=False
        )
        plt.xlabel("FDL1")
        plt.ylabel("FDL2")
        plt.title(f"{pt_to_bcpall[sample_name]}", fontsize=14)
        plt.savefig(os.path.join(figdir, f"{sample_name}_{suffix}.png"), dpi=300, bbox_inches="tight")
        plt.close()

    # ---------------------
    # Density histogram
    # ---------------------
    plt.figure()
    plt.hist(ad.obs["mellon_log_density"], bins=100)
    plt.xlabel("Log density")
    plt.ylabel("Cell count")
    plt.title(f"{pt_to_bcpall[sample_name]}", fontsize=14)
    plt.savefig(os.path.join(figdir, f"{sample_name}_density_histogram.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # ---------------------
    # Save processed AnnData
    # ---------------------
    out_h5ad = os.path.join(OUT_H5AD_DIR, f"{sample_name}_mellon_processed.h5ad")
    ad.write(out_h5ad)
    logger.info(f"Saved processed AnnData: {out_h5ad}")
    logger.info(f"=== Finished {sample_name} ===")

# =========================
# Entry point
# =========================

if __name__ == "__main__":
    logger.info("Starting Mellon batch analysis")
    for sample_name, cfg in SAMPLES.items():
        run_mellon_pipeline(sample_name, cfg)
    logger.info("All samples processed successfully.")
    print("All done! Figures saved in sample folders.")
