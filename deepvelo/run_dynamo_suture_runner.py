#!/usr/bin/env python3
"""
Run Dynamo analysis on scvelo/E17_suture_velocity.h5ad with robust handling.
"""
import argparse
import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

try:
    import dynamo as dyn
except Exception as exc:
    raise ImportError(
        "dynamo is required to run this script. "
        "Activate the correct environment and install dynamo."
    ) from exc

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = ROOT / "scvelo" / "E17_suture_velocity.h5ad"
DEFAULT_OUTPUT = ROOT / "dynamo" / "e17_suture"
DEFAULT_OUTPUT.mkdir(parents=True, exist_ok=True)

CLUSTER_COLORS = {
    "OG1": "#FF0000",
    "OG2": "#FF8000",
    "OG3": "#FFFF00",
    "OG4": "#00FF00",
    "PO1": "#0080FF",
    "PO2": "#8000FF",
}
DEFAULT_COLOR = "#AAAAAA"
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
            if any(str(v).startswith("OG") or str(v).startswith("PO") for v in vals):
                return col
    for col in adata.obs.columns:
        vals = adata.obs[col].astype(str).unique()
        if any(str(v).startswith("OG") or str(v).startswith("PO") for v in vals):
            return col
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if pd.api.types.is_categorical_dtype(dtype) or dtype == object:
            return col
    return None


def ensure_umap(adata, n_pcs=30):
    if "X_umap" in adata.obsm:
        return
    print("  - UMAP not found; computing PCA, neighbors, and UMAP.")
    if "X_pca" not in adata.obsm:
        if adata.shape[1] < 2:
            raise RuntimeError("Not enough genes to compute PCA.")
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(n_pcs, adata.shape[1] - 1))
    sc.pp.neighbors(adata, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata)


def plot_umap_by_cluster(adata, cluster_col, out_path):
    umap = adata.obsm.get("X_umap")
    if umap is None:
        raise RuntimeError("UMAP coordinates missing; cannot plot UMAP.")
    labels = (
        adata.obs[cluster_col].astype(str)
        if cluster_col
        else pd.Series(["unknown"] * adata.n_obs, index=adata.obs_names)
    )
    unique_labels = sorted(labels.unique(), key=str)
    fig, ax = plt.subplots(figsize=(10, 8))
    for label in unique_labels:
        mask = labels == label
        color = CLUSTER_COLORS.get(label, DEFAULT_COLOR)
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=25,
            c=color,
            alpha=0.7,
            label=label,
            linewidths=0,
        )
    ax.set_title("Dynamo UMAP by Cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), title="Cluster")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_velocity_vectors(adata):
    if "velocity_umap" in adata.obsm:
        vel = adata.obsm["velocity_umap"]
        if vel.shape[0] == adata.n_obs and vel.shape[1] >= 2:
            return vel[:, :2]
    for layer_name in ["velocity", "velocity_u", "velocity_umap"]:
        if layer_name in adata.layers:
            vel = adata.layers[layer_name]
            if vel.shape[0] == adata.n_obs and vel.shape[1] >= 2:
                return vel[:, :2]
    for entry in adata.obsm.keys():
        if "vel" in entry.lower() or "velocity" in entry.lower():
            vel = adata.obsm[entry]
            if vel.shape[0] == adata.n_obs and vel.shape[1] >= 2:
                return vel[:, :2]
    return None


def plot_velocity_quiver(adata, out_path, cluster_col=None, subset=500, scale=1.0):
    umap = adata.obsm.get("X_umap")
    if umap is None:
        raise RuntimeError("UMAP coordinates missing; cannot generate velocity quiver plot.")
    vel = get_velocity_vectors(adata)
    if vel is None:
        raise RuntimeError("No velocity vectors found for quiver plot.")
    n_cells = adata.n_obs
    idx = np.random.choice(n_cells, subset, replace=False) if n_cells > subset else np.arange(n_cells)
    fig, ax = plt.subplots(figsize=(10, 8))
    if cluster_col and cluster_col in adata.obs:
        labels = adata.obs[cluster_col].astype(str)
        for label in sorted(labels.unique(), key=str):
            mask = labels == label
            color = CLUSTER_COLORS.get(label, DEFAULT_COLOR)
            ax.scatter(
                umap[mask, 0],
                umap[mask, 1],
                s=15,
                c=color,
                alpha=0.25,
                linewidths=0,
            )
    else:
        ax.scatter(umap[:, 0], umap[:, 1], s=15, c="#CCCCCC", alpha=0.25, linewidths=0)
    ax.quiver(
        umap[idx, 0],
        umap[idx, 1],
        vel[idx, 0],
        vel[idx, 1],
        angles="xy",
        scale_units="xy",
        scale=1.0 / scale,
        width=0.003,
        color="#333333",
        alpha=0.8,
    )
    ax.set_title("Dynamo Velocity Quiver on UMAP")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def process_suture_file(input_path: Path, output_dir: Path, cores: int, cluster_column: str = None):
    print(f"Processing {input_path}")
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    adata = sc.read_h5ad(str(input_path))
    print(f"  - Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    if cluster_column and cluster_column not in adata.obs:
        print(f"  - Requested cluster column '{cluster_column}' not found; autodetecting instead.")
        cluster_column = None
    if cluster_column is None:
        cluster_column = find_cluster_column(adata)
    print(f"  - Cluster column: {cluster_column}")
    ensure_umap(adata)
    try:
        print("  - Running dynamo.pp.recipe_velocity()")
        dyn.pp.recipe_velocity(adata)
    except Exception:
        print("  - Warning: dynamo.pp.recipe_velocity() failed; continuing.")
        traceback.print_exc()
    try:
        print("  - Running dynamo.tl.dynamics()")
        dyn.tl.dynamics(adata, cores=cores, verbose=False)
    except Exception:
        print("  - Warning: dynamo.tl.dynamics() failed.")
        traceback.print_exc()
    try:
        print("  - Running dynamo.tl.cell_velocities()")
        dyn.tl.cell_velocities(adata, basis="umap", use_jacobian=False)
    except Exception:
        print("  - Warning: dynamo.tl.cell_velocities() failed.")
        traceback.print_exc()
    try:
        print("  - Building dynamo.vf.VectorField()")
        dyn.vf.VectorField(adata, basis="umap", cores=cores, verbose=False)
    except Exception:
        print("  - Warning: dynamo.vf.VectorField() failed.")
        traceback.print_exc()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = output_dir / f"{input_path.stem}_dynamo_processed.h5ad"
    print(f"  - Saving processed file to {out_h5ad}")
    adata.write_h5ad(str(out_h5ad))
    plot_umap_by_cluster(adata, cluster_column, output_dir / f"{input_path.stem}_umap_clusters.png")
    try:
        plot_velocity_quiver(adata, output_dir / f"{input_path.stem}_velocity_quiver.png", cluster_column)
    except Exception as exc:
        print(f"  - Warning: velocity quiver plot failed: {exc}")
    print(f"Finished processing {input_path.name}")


def parse_args():
    parser = argparse.ArgumentParser(description="Run Dynamo analysis on a scVelo h5ad file.")
    parser.add_argument(
        "--input-file",
        type=Path,
        default=DEFAULT_INPUT,
        help="Path to scvelo/E17_suture_velocity.h5ad",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Directory to write Dynamo results",
    )
    parser.add_argument(
        "--cluster-column",
        type=str,
        default=None,
        help="Cluster column name to use for UMAP coloring",
    )
    parser.add_argument("--cores", type=int, default=4, help="Number of cores for Dynamo computations")
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        process_suture_file(args.input_file, args.output_dir, args.cores, args.cluster_column)
        print(f"\nDynamo analysis finished. Outputs written to {args.output_dir}")
    except Exception:
        print("Fatal error during processing.")
        traceback.print_exc()


if __name__ == "__main__":
    main()