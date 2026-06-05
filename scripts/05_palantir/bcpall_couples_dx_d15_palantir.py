#!/usr/bin/env python
# coding: utf-8

"""
Palantir analysis on multiple BCP-ALL mmSW samples
- save-only plotting
- batch processing of multiple h5ad files
"""

# =========================
# Imports & global setup
# =========================

import palantir
import scanpy as sc
import pandas as pd
import os
import re

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import warnings
from numba.core.errors import NumbaDeprecationWarning

# ---- Silence warnings ----
warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

from tqdm import tqdm

# ---- Inline plotting ----
try:
    from IPython import get_ipython
    ip = get_ipython()
    if ip is not None:
        ip.run_line_magic('matplotlib', 'inline')
except Exception:
    # Not running inside IPython/Jupyter — skip inline backend setup
    pass

# ---- Matplotlib: save-only, no GUI ----
import logging
import matplotlib
matplotlib.use("Agg") 
import matplotlib.font_manager as font_manager

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
    "PT06": {
        "h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt06.h5ad",
        "start_cell": "r2_326029",
        "terminal_states": {
            "r2_299597": "Mono",
            "r2_465074": "B cell",
        },
    },
    "PT07": {
        "h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt07.h5ad",
        "start_cell": "r2_619078",
        "terminal_states": {
            "r2_12842918": "Mono",
            "r2_1779925": "B cell",
        },
    },
    "PT09": {
        "h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt09.h5ad",
        "start_cell": "r4_12992684",
        "terminal_states": {
            "r4_10345014": "Mono",
            "r4_461964": "B cell",
        },
    },
    "PT10": {
        "h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt10.h5ad",
        "start_cell": "r4_10645650",
        "terminal_states": {
            "r4_167061": "Mono",
            "r4_13155111": "B cell",
        },
    },
    "PT11": {
        "h5ad": "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/h5ad_objects_merge_dx/pt11.h5ad",
        "start_cell": "r5_3858494",
        "terminal_states": {
            "r5_33478": "Mono",
            "r5_162122": "B cell",
        },
    },
}

OUT_H5AD_DIR = (
    "/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/palantir/processed_palantir_h5ad")

FIGURES_BASE_DIR = ("/media/data/lab/BCPALL-MMSwitching-SingleCellLandscape/palantir/figures/")

GENES_TRENDS = ["RAG1", "MME", "AGAP1", "CD33", "CEBPD", "PLAGL1"]

# =========================
# Pipeline
# =========================

def run_palantir_pipeline(sample_name: str, cfg: dict):

    h5ad_path = cfg["h5ad"]
    start_cell = cfg["start_cell"]
    terminal_states = pd.Series(cfg["terminal_states"])

    # Output dirs
    figdir = os.path.join(FIGURES_BASE_DIR, sample_name)
    os.makedirs(figdir, exist_ok=True)
    sc.settings.figdir = figdir
    os.makedirs(OUT_H5AD_DIR, exist_ok=True)

    # Loading data

    logger.info(f"=== Running Palantir for {sample_name} ===")

    ad = sc.read(h5ad_path)
    ad

    # Run diffusion maps
    dm_res = palantir.utils.run_diffusion_maps(ad, pca_key="X_pca", n_components=5)

    ms_data = palantir.utils.determine_multiscale_space(ad)


    # Visualization

    sc.pp.neighbors(ad)
    sc.tl.draw_graph(ad)
    sc.set_figure_params(dpi_save=300)


    # Use scanpy functions to visualize umaps or FDL

    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        show=False,
        frameon=False,
    )
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_draw_graph_fa.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    sample_name_colors = {
        "BCPALL#06_DX" : "#d992eb",
        "BCPALL#06_D15" : "#AC92EB",
        "BCPALL#07_DX" : "#4FC1E8",
        "BCPALL#07_D15" : "#4f75e8",
        "BCPALL#08_DX" : "#94b76c",
        "BCPALL#08_D15" : "#607f3e",
        "BCPALL#09_DX" : "#FFCE54",
        "BCPALL#09_D15" : "#d39700",
        "BCPALL#10_DX" : "#ED5564",
        "BCPALL#10_D15" : "#c71627",
        "BCPALL#11_DX" : "#ED9255",
        "BCPALL#11_D15" : "#985726"
    }

    # Create mapping from old names to new names
    name_mapping = {
        "PT06_DX": "BCPALL#06_DX",
        "PT06_D15": "BCPALL#06_D15",
        "PT07_DX": "BCPALL#07_DX",
        "PT07_D15": "BCPALL#07_D15",
        "PT08_DX": "BCPALL#08_DX",
        "PT08_D15": "BCPALL#08_D15",
        "PT13_DX": "BCPALL#09_DX",
        "PT13_D15": "BCPALL#09_D15",
        "PT14_DX": "BCPALL#10_DX",
        "PT14_D15": "BCPALL#10_D15",
        "PT15_DX": "BCPALL#11_DX",
        "PT15_D15": "BCPALL#11_D15",
    }

    # Create reverse mapping for patient names
    pt_to_bcpall = {
        "PT06": "BCPALL#06",
        "PT07": "BCPALL#07",
        "PT08": "BCPALL#08",
        "PT13": "BCPALL#09",
        "PT14": "BCPALL#10",
        "PT15": "BCPALL#11",
    }

    # Rename categories in observations
    ad.obs['Sample_Name'] = ad.obs['Sample_Name'].cat.rename_categories(name_mapping)

    # Filter colors to only show the current patient
    patient_bcpall = pt_to_bcpall[sample_name]
    filtered_categories = [cat for cat in ad.obs['Sample_Name'].cat.categories if cat.startswith(patient_bcpall)]
    ad.obs['Sample_Name'] = ad.obs['Sample_Name'].cat.set_categories(filtered_categories)
    ad.uns["Sample_Name_colors"] = [sample_name_colors[cat] for cat in filtered_categories]


    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "normal"
    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color="Sample_Name",
        frameon=True,
        title='',
        wspace=4,
        legend_fontweight='normal',
        legend_fontoutline=20,
        show=False,
    )
    plt.xlabel("FDL1")
    plt.ylabel("FDL2")
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_color_sample.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


    # MAGIC imputation 

    imputed_X = palantir.utils.run_magic_imputation(ad)

    # add plotting settings and graphics

    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "italic"
    plt.figure(figsize=(10,8))
    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color=["CD34", "CD38","CD33",  "CLEC12A", "CD2", "CD58", "PTPRC","CD19","MS4A1","CLEC4E", "CD24","MME","CEBPA", "CEBPD","CD14","FCGR3A","NR3C1", "CD22", "CD44", "FLT3", "TGFB1", "EBF1", 
            "CD79B", "RAG1", "IFNB1", "IRF8"],
        frameon=False,
        show=False,
    )
    plt.xlabel("FDL1")
    plt.ylabel("FDL2")
    plt.tight_layout()
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_genes_markers.png"),
        dpi=300,
    )
    plt.close()

    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color=["dux4_up_score1"],
        show=False,
        frameon=False,
    )

    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "normal"
    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color=["dux4_up_score1"],
        frameon=True,
        wspace=4,
        title=f'{pt_to_bcpall[sample_name]}',
        legend_fontweight='normal',
        legend_fontoutline=20,
        show=False,
    )
    plt.xlabel("FDL1")
    plt.ylabel("FDL2")
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_dux4_score.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    

    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "normal"
    plot_genes_12 = [
        "PAX5", "CD34", "AGAP1", "EBF1", "MME", "PIK3R1",
        "IGLL1", "SOX4", "PLAGL1", "VPREB1", "TGFB1", "CD47"
    ]
    for gene in plot_genes_12:
        plt.figure(figsize=(6,5))
        sc.pl.embedding(
            ad,
            basis="X_draw_graph_fa",
            layer="MAGIC_imputed_data",
            color=[gene],
            title="",
            frameon=True,
            wspace=4,
            legend_fontweight='normal',
            legend_fontoutline=20,
            show=False,
        )
        plt.xlabel("FDL1")
        plt.ylabel("FDL2")
        ax = plt.gca()
        ax.set_title(gene, fontstyle='italic', fontfamily='Karla', fontsize=24)
        plt.tight_layout()
        plt.savefig(
            os.path.join(figdir, f"{sample_name}_{gene}.png"),
            dpi=300,
        )
        plt.close()

    # ---- Additional: plot NEW_GENES into a separate folder per sample ----
    if NEW_GENES:
        new_dir = os.path.join(NEW_GENES_BASE_DIR, sample_name)
        os.makedirs(new_dir, exist_ok=True)
        for gene in NEW_GENES:
            plt.figure(figsize=(6,5))
            sc.pl.embedding(
                ad,
                basis="X_draw_graph_fa",
                layer="MAGIC_imputed_data",
                color=[gene],
                title="",
                frameon=True,
                show=False,
            )
            ax = plt.gca()
            ax.set_title(gene, fontstyle='italic', fontfamily='Karla', fontsize=18)
            ax.set_xlabel("FDL1")
            ax.set_ylabel("FDL2")
            plt.tight_layout()
            plt.savefig(
                os.path.join(new_dir, f"{sample_name}_{gene}.png"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

    # --- Ensure a consistent color mapping for predicted_CellType_BoneMarrowMap ---
    # Use the full list of known bone-marrow celltype labels (collected from
    # the dataset) and assign a deterministic palette so the legend and colors
    # are identical across patients. If you want a different ordering, edit
    # DEFAULT_BM_CELLTYPES below.
    DEFAULT_BM_CELLTYPES = [
        'ASDC', 'BFU-E', 'Basophilic Erythroblast', 'CD14 Mono', 'CD16 Mono',
        'CD4 Central Memory', 'CD4 Effector Memory', 'CD4 Naive',
        'CD8 Central Memory', 'CD8 Effector Memory 1', 'CD8 Effector Memory 2',
        'CD8 Tissue Resident Memory', 'CLP', 'Cycling Progenitor', 'Early GMP',
        'Early ProMono', 'EoBasoMast Precursor', 'GMP-Mono', 'HSC', 'Immature B',
        'LMPP', 'Large Pre-B', 'Late ProMono', 'MEP', 'MLP', 'MLP-II',
        'MPP-MkEry', 'MPP-MyLy', 'Mature B', 'Megakaryocyte Precursor',
        'NK', 'NK CD56high', 'Orthochromatic Erythroblast',
        'Polychromatic Erythroblast', 'Pre-ProB', 'Pre-cDC', 'Pre-pDC',
        'Pre-pDC Cycling', 'Pro-B Cycling', 'Pro-B VDJ', 'Small Pre-B',
        'Stromal', 'T Proliferating', 'cDC1', 'cDC2', 'pDC'
    ]

    # Build the final category list (keep DEFAULT order, then any extras found
    # in the current AnnData that are not in the default list).
    final_cats = list(DEFAULT_BM_CELLTYPES)
    if "predicted_CellType_BoneMarrowMap" in ad.obs.columns:
        # Treat explicit 'NA' string as missing and remove it from categories
        ad.obs["predicted_CellType_BoneMarrowMap"] = ad.obs["predicted_CellType_BoneMarrowMap"].replace('NA', pd.NA)
        present = [str(x) for x in ad.obs["predicted_CellType_BoneMarrowMap"].dropna().unique()]
        # filter out any literal 'NA' entries if present
        present = [p for p in present if p != 'NA']
        extras = [c for c in sorted(present) if c not in final_cats]
        final_cats.extend(extras)

    # Generate a deterministic palette with N colors (one per category).
    # Use a continuous categorical cmap so it scales for >20 categories.
    import matplotlib as mpl
    cmap = mpl.cm.get_cmap("tab20", len(final_cats))
    celltype_colors = {cat: mpl.colors.to_hex(cmap(i)) for i, cat in enumerate(final_cats)}

    # Make the obs column categorical with the full set of categories so the
    # legend shows the same entries (absent categories will simply have no
    # points plotted for this sample).
    if "predicted_CellType_BoneMarrowMap" in ad.obs.columns:
        ad.obs["predicted_CellType_BoneMarrowMap"] = (
            ad.obs["predicted_CellType_BoneMarrowMap"].astype("category").cat.set_categories(final_cats)
        )

    # Tell Scanpy which colors to use (in the same order as categories)
    # Remove any 'NA' from the final mapping just in case, and build colors list
    final_cats = [c for c in final_cats if c != 'NA']
    ad.uns["predicted_CellType_BoneMarrowMap_colors"] = [celltype_colors.get(cat, "#BBBBBB") for cat in final_cats]

    # Plot using the fixed mapping
    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color="predicted_CellType_BoneMarrowMap",
        show=False,
        frameon=True,
    )
    # Ensure axes are visible (Scanpy sometimes hides them). Find the Axes
    # that contains the scatter (PathCollection) and enable its axis/labels.
    fig = plt.gcf()
    ax = None
    for candidate in fig.get_axes():
        # pick the axes that actually has plotted collections (the embedding)
        if any(hasattr(col, "get_offsets") for col in candidate.collections):
            ax = candidate
            break
    if ax is None:
        # fallback to first axis or current axis
        axes = fig.get_axes()
        ax = axes[0] if axes else plt.gca()

    ax.axis("on")
    ax.set_xlabel("FDL1")
    ax.set_ylabel("FDL2")
    plt.tight_layout()
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_predicted_CellType_BoneMarrowMap.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    

    sc.pl.embedding(
        ad,
        basis="X_draw_graph_fa",
        layer="MAGIC_imputed_data",
        color="seurat_clusters",
        show=False,
        frameon=False,
    )
    plt.savefig(
        os.path.join(figdir, f"{sample_name}_seurat_clusters.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    


    # Running Palantir

    terminal_states = pd.Series(cfg["terminal_states"])
    start_cell_id = start_cell
    sc.pp.neighbors(ad, use_rep="X_pca")

    palantir.utils.run_diffusion_maps(
        ad,
        pca_key="X_pca",
        n_components=5,
    )

    palantir.utils.determine_multiscale_space(ad)

    start_cell_plot = pd.Series(
        [start_cell_id],
        index=[start_cell_id],
    )

    palantir.plot.highlight_cells_on_umap(ad, terminal_states, embedding_basis='X_draw_graph_fa')
    palantir.plot.highlight_cells_on_umap(ad, start_cell_plot, embedding_basis='X_draw_graph_fa')
   
    pr_res = palantir.core.run_palantir(
        ad, start_cell_id, num_waypoints=500, terminal_states=terminal_states
    )

    # Visualizing Palantir results
    print(ad.obs.columns)

    # Create temporary PResults without branch probabilities
    pr_plot = palantir.presults.PResults(
        pr_res.pseudotime,
        pr_res.entropy,
        pd.DataFrame(index=pr_res.branch_probs.index),
        None,
    )
    palantir.plot.plot_palantir_results(
        ad,
        pr_res=pr_plot,
        s=3,
        embedding_basis="X_draw_graph_fa",
    )

    # --- Axis formatting ---
    fig = plt.gcf()
    axes = fig.get_axes()

    # Scatter plot axes (axes[0]=pseudotime, axes[2]=entropy)
    for ax in [axes[0], axes[2]]:
        ax.axis("on")
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("FDL1")
        ax.set_ylabel("FDL2")

    # Colorbar axes (axes[1]=pseudotime colorbar, axes[3]=entropy colorbar)
    for ax in [axes[1], axes[3]]:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.yaxis.set_major_locator(ticker.AutoLocator())
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.tick_params(axis="y", which="both", left=False, labelleft=False,
                       right=True, labelright=True, length=3)

    plt.savefig(
        os.path.join(figdir, f"palantir_{sample_name}.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    # Keep only selected fate probabilities
    ad.obsm["palantir_fate_probabilities"] = (
        ad.obsm["palantir_fate_probabilities"][["Mono", "B cell"]]
    )
    print(ad.obsm["palantir_fate_probabilities"])

    # ---- Trajectories ----
    for fate in terminal_states.values:
        masks = palantir.presults.select_branch_cells(ad, q=.01, eps=.01)
        palantir.plot.plot_trajectory(
            ad,
            fate,
            embedding_basis="X_draw_graph_fa",
            smoothness=0.002,
        )
        plt.savefig(
            os.path.join(figdir, f"{fate}_trajectory_{sample_name}.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    # Visualizing the branch selection

    palantir.plot.plot_branch_selection(ad, embedding_basis='X_draw_graph_fa')

    palantir.plot.plot_trajectory(ad, "Mono", embedding_basis='X_draw_graph_fa', smoothness=0.002)
    plt.savefig(
        os.path.join(figdir, f"mono_branch_{sample_name}.png"),
        dpi=300,
        bbox_inches="tight",
    )    
    plt.close()

    palantir.plot.plot_trajectory(ad, "B cell", embedding_basis='X_draw_graph_fa', smoothness=0.002)
    plt.savefig(
        os.path.join(figdir, f"bcell_branch_{sample_name}.png"),
        dpi=300,
        bbox_inches="tight",
    )  
    plt.close()

    # ---- Trajectories with consistent figure size, font, and labels ----
    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "normal"
    plt.rcParams["font.size"] = 14          
    plt.rcParams["axes.labelsize"] = 18 
    plt.rcParams["axes.titlesize"] = 20     

    for fate in terminal_states.values:
        plt.figure(figsize=(6,5))
        palantir.plot.plot_trajectory(
            ad,
            fate,
            embedding_basis="X_draw_graph_fa",
            smoothness=0.002,
        )
        fig = plt.gcf()
        ax = fig.get_axes()[0]
        ax.axis("on")
        ax.set_xlabel("FDL1")
        ax.set_ylabel("FDL2")
        plt.tight_layout()
        plt.savefig(
            os.path.join(figdir, f"{fate}_trajectory_{sample_name}.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    # Mono and B cell branches
    for fate in ["Mono", "B cell"]:
        plt.figure(figsize=(6,5))
        palantir.plot.plot_trajectory(
            ad,
            fate,
            embedding_basis='X_draw_graph_fa',
            smoothness=0.002
        )
        fig = plt.gcf()
        ax = fig.get_axes()[0]
        ax.axis("on")
        ax.set_xlabel("FDL1")
        ax.set_ylabel("FDL2")
        plt.tight_layout()
        plt.savefig(
            os.path.join(figdir, f"{fate.lower()}_branch_{sample_name}.png"),
            dpi=300,
            bbox_inches="tight",
        )    
        plt.close()
    
    # Gene trends in 2x3 grid
    palantir.presults.compute_gene_trends(
        ad,
        expression_key="MAGIC_imputed_data",
    )

    plt.rcParams["font.family"] = "Karla"
    plt.rcParams["font.style"] = "normal"
    plt.rcParams["font.size"] = 14

        # ---- Choose grid layout depending on sample ----
    if sample_name == "PT13":  # BCPALL#09
        nrows, ncols = 2, 3
        figsize = (16, 6)
    else:
        nrows, ncols = 3, 2
        figsize = (10, 10)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten()

    legend_handles, legend_labels = [], []

    for ax, gene in zip(axes, GENES_TRENDS):
        gene_fig = palantir.plot.plot_gene_trends(ad, [gene])
        gene_ax = gene_fig.get_axes()[0]

        for line in gene_ax.get_lines():
            l, = ax.plot(line.get_xdata(), line.get_ydata(),
                    color=line.get_color(), lw=line.get_linewidth(), label=line.get_label())
            if line.get_label() not in legend_labels:
                legend_handles.append(l)
                legend_labels.append(line.get_label())
        for collection in gene_ax.collections:
            ax.add_collection(collection)
        ax.autoscale()

        ax.set_title(gene, fontsize=18, fontstyle="italic", fontfamily="Karla")
        ax.set_xlabel("Pseudotime", fontsize=14, fontfamily="Karla")
        ax.set_ylabel("Expression", fontsize=14, fontfamily="Karla")
        ax.grid(False)
        ax.legend().remove()
        plt.close(gene_fig)

    leg=fig.legend(
        legend_handles, legend_labels,
        loc="lower center",
        ncol=len(legend_labels),
        fontsize=16,
        frameon=False,
        bbox_to_anchor=(0.5, -0.04),
        handlelength=1.5,
    )
    for legline in leg.get_lines():
        legline.set_linewidth(2)

    plt.tight_layout()
    plt.savefig(
        os.path.join(figdir, f"genes_trends_{sample_name}.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    # Save processed AnnData

    out_h5ad = os.path.join(
        OUT_H5AD_DIR, f"{sample_name}_palantir_processed.h5ad"
    )
    ad.write(out_h5ad)

    logger.info(f"=== Finished {sample_name} ===")

# =========================
# Entry point
# =========================

if __name__ == "__main__":
    logger.info("Starting Palantir batch analysis")

    for sample_name, cfg in SAMPLES.items():
        run_palantir_pipeline(sample_name, cfg)

    logger.info("All samples processed successfully.")

    print("All done! Figures are saved in folders ./figures. Bye!...")
