"""
修正版：ストリームプロット + CellRank
APIの問題に対応
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

import scanpy as sc
import scvelo as scv
import cellrank as cr

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

print("Loading data...")
adata = sc.read_h5ad("E17_og_integrated_harmony_filtered_velocity.h5ad")
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

outdir = "figures_integrated_analysis"
os.makedirs(outdir, exist_ok=True)
sc.settings.figdir = outdir

OG_COLORS = {
    "OG1": "#e41a1c",
    "OG2": "#ff7f00",
    "OG3": "#4daf4a",
    "OG4": "#984ea3",
    "PO1": "#377eb8",
    "PO2": "#a65628",
}

# ============================================================
# 1. Stream plot (alternative method)
# ============================================================
print("\n=== Generating Stream Plot ===")
try:
    # Method: use velocity_embedding_grid
    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="og_cluster",
        palette=OG_COLORS,
        title="RNA Velocity (Grid, Batch Corrected)",
        save="_velocity_grid.png",
        dpi=150,
    )
    print("Grid plot saved.")
except Exception as e:
    print(f"Grid plot error: {e}")

try:
    # Method: just scatter with velocity arrows
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.embedding(adata, basis="umap", color="og_cluster", palette=OG_COLORS, ax=ax, show=False)
    # Add quiver plot manually
    if "velocity_umap" in adata.obsm:
        X = adata.obsm["umap"]
        V = adata.obsm["velocity_umap"]
        ax.quiver(X[:, 0], X[:, 1], V[:, 0], V[:, 1], alpha=0.3, scale=10, width=0.002)
    ax.set_title("RNA Velocity (Manual Quiver)")
    plt.tight_layout()
    plt.savefig(f"{outdir}/velocity_quiver.png", dpi=150)
    plt.close()
    print("Quiver plot saved.")
except Exception as e:
    print(f"Quiver plot error: {e}")

# ============================================================
# 2. CellRank with corrected API
# ============================================================
print("\n=== CellRank Analysis ===")
try:
    # Build kernels
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    
    # Combine kernels properly using correct method
    combined = 0.8 * vk + 0.2 * ck
    
    # Use cr.tl functions instead of kernel methods
    print("Computing macrostates...")
    cr.tl.macrostates(combined, n_states=4)
    print(f"Macrostates computed: {adata.obs['macrostates'].value_counts().to_dict()}")
    
    print("Computing terminal states...")
    cr.tl.terminal_states(combined)
    print(f"Terminal states: {adata.obs['terminal_states'].value_counts().to_dict()}")
    
    print("Computing absorption probabilities...")
    cr.tl.absorption_probabilities(combined)
    print("Absorption probabilities computed.")
    
    print("CellRank analysis completed.")
    
    # Visualization
    try:
        sc.pl.umap(adata, color="macrostates", palette="tab20", save="_macrostates.png")
        print("Macrostates plot saved.")
    except:
        pass
    
    try:
        sc.pl.umap(adata, color="terminal_states", palette="tab20", save="_terminal_states.png")
        print("Terminal states plot saved.")
    except:
        pass
    
    # Print results
    if "absorption_probabilities" in adata.obsm:
        ab_probs = adata.obsm["absorption_probabilities"]
        print(f"Absorption probabilities shape: {ab_probs.shape}")
        print(f"Terminal states: {ab_probs.index.tolist()}")
    
except Exception as e:
    print(f"CellRank error: {e}")
    import traceback
    traceback.print_exc()

# ============================================================
# 3. Save results
# ============================================================
adata.write(f"{outdir}/analysis_results_fixed.h5ad")
print(f"\nAnalysis completed. Results saved to {outdir}/")
