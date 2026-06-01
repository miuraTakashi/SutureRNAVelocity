"""
ストリームプロット（streamlines）の生成
"""

import os
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import griddata

warnings.filterwarnings("ignore")

import scanpy as sc
import scvelo as scv

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

print("Loading data...")
adata = sc.read_h5ad("E17_og_integrated_harmony_filtered_velocity.h5ad")

outdir = "figures_integrated_analysis"
os.makedirs(outdir, exist_ok=True)

OG_COLORS = {
    "OG1": "#e41a1c",
    "OG2": "#ff7f00",
    "OG3": "#4daf4a",
    "OG4": "#984ea3",
    "PO1": "#377eb8",
    "PO2": "#a65628",
}

# ============================================================
# Method 1: Use scVelo's built-in streamline function
# ============================================================
print("\n=== Attempting scVelo streamlines ===")
try:
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="og_cluster",
        palette=OG_COLORS,
        title="RNA Velocity Stream (Batch Corrected OG/PO)",
        density=1.5,
        linewidth=0.5,
        show=False,
        save="_stream.png",
    )
    print("scVelo streamlines succeeded!")
except Exception as e:
    print(f"scVelo streamlines failed: {e}")
    
    # ============================================================
    # Method 2: Manual streamline implementation
    # ============================================================
    print("\n=== Using manual streamline implementation ===")
    try:
        umap_coords = adata.obsm["X_umap"]
        velocity = adata.layers["velocity"]
        
        # Compute velocity in UMAP space
        # Use velocity + PCA projection
        if "velocity_umap" not in adata.obsm:
            print("Computing velocity_umap...")
            scv.tl.velocity_embedding(adata, basis="umap")
        
        V_umap = adata.obsm["velocity_umap"]
        
        # Create grid for streamplot
        x = np.linspace(umap_coords[:, 0].min() - 0.5, umap_coords[:, 0].max() + 0.5, 30)
        y = np.linspace(umap_coords[:, 1].min() - 0.5, umap_coords[:, 1].max() + 0.5, 30)
        X, Y = np.meshgrid(x, y)
        
        # Interpolate velocity vectors on grid
        U = griddata(umap_coords, V_umap[:, 0], (X, Y), method='cubic', fill_value=0)
        V = griddata(umap_coords, V_umap[:, 1], (X, Y), method='cubic', fill_value=0)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot cells colored by cluster
        colors_array = [OG_COLORS.get(c, "#cccccc") for c in adata.obs["og_cluster"]]
        ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors_array, s=20, alpha=0.6, edgecolors="none")
        
        # Plot streamlines
        speed = np.sqrt(U**2 + V**2)
        lw = 2 * speed / speed.max()  # Normalize line width
        
        strm = ax.streamplot(X, Y, U, V, linewidth=lw, color="black", density=1.5, arrowsize=1.5)
        
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title("RNA Velocity Stream (Batch Corrected OG/PO)")
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=OG_COLORS[c], label=c) for c in sorted(OG_COLORS.keys())]
        ax.legend(handles=legend_elements, loc="best", title="Cluster", fontsize=9)
        
        plt.tight_layout()
        plt.savefig(f"{outdir}/scvelo__stream.png", dpi=150, bbox_inches="tight")
        plt.close()
        
        print("Manual streamlines saved!")
        
    except Exception as e2:
        print(f"Manual streamlines failed: {e2}")
        import traceback
        traceback.print_exc()

print("\nDone.")
