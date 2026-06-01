"""
最終実装：簡潔な分析
- サイクル検出（既に完了）
- 基本的な可視化
- 結果の保存
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

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

print("Loading data...")
adata = sc.read_h5ad("E17_og_integrated_harmony_filtered_velocity.h5ad")
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

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
# 1. UMAP + Velocity overlay
# ============================================================
print("\n=== Creating UMAP with Velocity Overlay ===")
try:
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # UMAP scatter
    umap_coords = adata.obsm.get("X_umap", adata.obsm.get("umap", None))
    if umap_coords is None:
        print("Warning: UMAP not found. Computing...")
        sc.tl.umap(adata)
        umap_coords = adata.obsm["X_umap"]
    
    colors = [OG_COLORS.get(c, "#cccccc") for c in adata.obs["og_cluster"]]
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=colors, s=30, alpha=0.7, edgecolors="none")
    
    # Velocity vectors
    if "velocity_umap" in adata.obsm:
        V = adata.obsm["velocity_umap"]
        ax.quiver(umap_coords[:, 0], umap_coords[:, 1], 
                 V[:, 0], V[:, 1], 
                 alpha=0.5, scale=20, width=0.003, headwidth=3, headlength=4)
    
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("RNA Velocity on UMAP (Batch Corrected OG/PO)")
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=OG_COLORS[c], label=c) for c in sorted(OG_COLORS.keys())]
    ax.legend(handles=legend_elements, loc="best", title="Cluster")
    
    plt.tight_layout()
    plt.savefig(f"{outdir}/velocity_overlay.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Velocity overlay saved.")
except Exception as e:
    print(f"Error: {e}")

# ============================================================
# 2. Cluster composition
# ============================================================
print("\n=== Cluster Composition ===")
composition = adata.obs["og_cluster"].value_counts().sort_index()
print(composition)

fig, ax = plt.subplots(figsize=(8, 6))
composition.plot(kind="bar", ax=ax, color=[OG_COLORS.get(c, "#cccccc") for c in composition.index])
ax.set_xlabel("Cluster")
ax.set_ylabel("Number of Cells")
ax.set_title("Cluster Composition")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{outdir}/cluster_composition.png", dpi=150)
plt.close()
print("Cluster composition saved.")

# ============================================================
# 3. Batch composition
# ============================================================
print("\n=== Batch Distribution ===")
batch_cluster = pd.crosstab(adata.obs["batch"], adata.obs["og_cluster"])
print(batch_cluster)

fig, ax = plt.subplots(figsize=(10, 6))
batch_cluster.plot(kind="bar", ax=ax, stacked=True)
ax.set_xlabel("Batch")
ax.set_ylabel("Number of Cells")
ax.set_title("Batch Distribution by Cluster")
plt.xticks(rotation=45)
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig(f"{outdir}/batch_distribution.png", dpi=150, bbox_inches="tight")
plt.close()
print("Batch distribution saved.")

# ============================================================
# 4. Summary statistics
# ============================================================
print("\n=== Summary Statistics ===")
summary = pd.DataFrame({
    "Total Cells": [adata.n_obs],
    "Total Genes": [adata.n_vars],
    "Clusters": [adata.obs["og_cluster"].nunique()],
    "Batches": [adata.obs["batch"].nunique()],
    "Has Velocity": ["velocity" in adata.layers],
    "Has Velocity Graph": ["velocity_graph" in adata.uns],
    "Has Velocity UMAP": ["velocity_umap" in adata.obsm],
})
print(summary.to_string(index=False))

# ============================================================
# 5. Save processed data and summary
# ============================================================
adata.write(f"{outdir}/analysis_results.h5ad")
summary.to_csv(f"{outdir}/summary_statistics.csv", index=False)
batch_cluster.to_csv(f"{outdir}/batch_cluster_distribution.csv")
composition.to_csv(f"{outdir}/cluster_composition.csv")

print(f"\nAll analyses completed. Results saved to {outdir}/")
print(f"  - velocity_overlay.png: UMAP with velocity vectors")
print(f"  - cluster_composition.png: Cell count per cluster")
print(f"  - batch_distribution.png: Batch composition")
print(f"  - analysis_results.h5ad: Full AnnData object")
print(f"  - summary_statistics.csv: Summary metrics")
