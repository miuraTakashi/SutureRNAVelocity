#!/usr/bin/env python3
"""
run_dynamo_analysis.py

Loads OG-filtered h5ad files from scvelo/, runs dynamo analysis (recipe_velocity,
dynamics, cell_velocities, VectorField), and saves plots and processed h5ad files
to dynamo/ folder.
"""
import os
import sys
import traceback
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import dynamo as dyn
except Exception as e:
    print("ERROR: Could not import dynamo. Install dynamo and try again.")
    raise

# Configuration
ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "scvelo"
OUTPUT_DIR = ROOT / "dynamo"
OUTPUT_DIR.mkdir(exist_ok=True)

BATCHES = ["E17_batch_1", "E17_batch_2", "E17_batch_3"]
SUFFIX = "_velocity_OG_filtered_v2.h5ad"

# cluster colors mapping
CLUSTER_COLORS = {
    "OG1": "#FF0000",
    "OG2": "#FF8000",
    "OG3": "#FFFF00",
    "OG4": "#00FF00",
    "PO1": "#0080FF",
    "PO2": "#8000FF",
}
DEFAULT_COLOR = "#AAAAAA"

# Candidate obs columns that may contain cluster labels
CANDIDATE_CLUSTER_COLS = [
    "cluster",
    "clusters",
    "louvain",
    "leiden",
    "seurat_clusters",
    "og_cluster",
    "OG_cluster",
    "assigned_cluster",
    "annotation",
    "group",
]


def find_cluster_column(adata):
    for col in CANDIDATE_CLUSTER_COLS:
        if col in adata.obs:
            vals = adata.obs[col].astype(str).unique()
            if any(v.startswith("OG") or v.startswith("PO") for v in vals):
                return col
    for col in adata.obs.columns:
        vals = adata.obs[col].astype(str).unique()
        if any(v.startswith("OG") or v.startswith("PO") for v in vals):
            return col
    for col in adata.obs.columns:
        if pd.api.types.is_categorical_dtype(adata.obs[col]) or adata.obs[col].dtype == object:
            return col
    return None


def ensure_umap(adata):
    if "X_umap" in adata.obsm:
        return
    print("  - UMAP not found; computing neighbors and UMAP.")
    if "X_pca" not in adata.obsm:
        print("  - PCA not found; computing PCA.")
        sc.pp.pca(adata)
    sc.pp.neighbors(adata, use_rep="X_pca" if "X_pca" in adata.obsm else None)
    sc.tl.umap(adata)


def plot_umap_by_cluster(adata, cluster_col, out_path, cluster_colors):
    umap = adata.obsm.get("X_umap")
    if umap is None:
        raise RuntimeError("UMAP coordinates missing in adata.obsm['X_umap']")
    df = pd.DataFrame(umap[:, :2], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    df["cluster"] = adata.obs[cluster_col].astype(str) if cluster_col else "unknown"
    unique = df["cluster"].unique()
    colors = []
    for u in unique:
        colors.append(cluster_colors.get(u, DEFAULT_COLOR))
    fig, ax = plt.subplots(figsize=(10, 8))
    for u, c in zip(unique, colors):
        sub = df[df["cluster"] == u]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=20, c=c, label=u, alpha=0.7, linewidths=0)
    ax.set_title(f"Dynamo UMAP by Cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fallback_quiver_from_velocity(adata, out_path, subset=500, scale=1.0):
    umap = adata.obsm.get("X_umap")
    if umap is None:
        raise RuntimeError("UMAP missing; cannot create quiver plot fallback.")
    velocity_keys = [k for k in adata.obsm.keys() if "vel" in k.lower() or "velocity" in k.lower()]
    vel = None
    for k in velocity_keys:
        v = adata.obsm.get(k)
        if v is None:
            continue
        if v.shape[0] == umap.shape[0] and v.shape[1] >= 2:
            vel = v[:, :2]
            break
    if vel is None:
        print("  - No velocity vector found in adata.obsm/layers for quiver fallback.")
        return False

    n_cells = umap.shape[0]
    if n_cells > subset:
        idx = np.random.choice(n_cells, subset, replace=False)
    else:
        idx = np.arange(n_cells)
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(umap[:, 0], umap[:, 1], s=10, c="#CCCCCC", alpha=0.5)
    Q = ax.quiver(umap[idx, 0], umap[idx, 1], vel[idx, 0], vel[idx, 1], angles="xy", scale_units="xy", scale=1.0/scale, width=0.003, color="#333333")
    ax.set_title("Dynamo Velocity Field (UMAP)")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def process_file(input_path: Path):
    print(f"\nProcessing {input_path.name}")
    try:
        adata = sc.read_h5ad(str(input_path))
    except Exception:
        print("  - ERROR: failed to read file.")
        traceback.print_exc()
        return

    try:
        print("  - Ensuring UMAP...")
        ensure_umap(adata)
    except Exception:
        print("  - ERROR: failed to compute/check UMAP.")
        traceback.print_exc()

    cluster_col = find_cluster_column(adata)
    print(f"  - Using cluster column: {cluster_col}")

    try:
        print("  - Running dynamo.pp.recipe_velocity()")
        dyn.pp.recipe_velocity(adata)
    except Exception:
        print("  - Warning: recipe_velocity failed; continuing anyway.")
        traceback.print_exc()

    try:
        print("  - Running dynamo.tl.dynamics()")
        dyn.tl.dynamics(adata, cores=4, verbose=False)
    except Exception:
        print("  - Warning: dynamo.tl.dynamics() failed.")
        traceback.print_exc()

    try:
        print("  - Running dynamo.tl.cell_velocities()")
        dyn.tl.cell_velocities(adata, basis="umap", use_jacobian=False)
    except Exception:
        print("  - Warning: dynamo.tl.cell_velocities() failed.")
        traceback.print_exc()

    vf_obj = None
    try:
        print("  - Building dynamo.vf.VectorField()")
        vf_obj = dyn.vf.VectorField(adata, basis="umap", cores=4, verbose=False)
    except Exception:
        print("  - Warning: VectorField construction failed.")
        traceback.print_exc()

    try:
        out_umap_cluster = OUTPUT_DIR / f"{input_path.stem}_umap_by_cluster.png"
        print(f"  - Plotting UMAP by cluster")
        plot_umap_by_cluster(adata, cluster_col, out_umap_cluster, CLUSTER_COLORS)
    except Exception:
        print("  - ERROR: failed to create UMAP cluster plot.")
        traceback.print_exc()

    try:
        out_quiver = OUTPUT_DIR / f"{input_path.stem}_velocity_quiver.png"
        ok = fallback_quiver_from_velocity(adata, out_quiver)
        if ok:
            print(f"  - Saved velocity quiver plot")
    except Exception:
        print("  - ERROR: failed to generate velocity quiver plot.")
        traceback.print_exc()

    try:
        out_h5ad = OUTPUT_DIR / f"{input_path.stem}_dynamo_processed.h5ad"
        print(f"  - Saving processed AnnData")
        adata.write_h5ad(str(out_h5ad))
    except Exception:
        print("  - ERROR: failed to save processed h5ad.")
        traceback.print_exc()

    print(f"Finished processing {input_path.name}")


def main():
    print("="*70)
    print("Starting dynamo batch analysis")
    print("="*70)
    for b in BATCHES:
        in_file = INPUT_DIR / f"{b}{SUFFIX}"
        if not in_file.exists():
            print(f"Skipping {b}: input file not found at {in_file}")
            continue
        try:
            process_file(in_file)
        except Exception:
            print(f"Unhandled error while processing {in_file.name}:")
            traceback.print_exc()

    print("\n" + "="*70)
    print("All done. Results are in:", OUTPUT_DIR)
    print("="*70)


if __name__ == "__main__":
    main()
