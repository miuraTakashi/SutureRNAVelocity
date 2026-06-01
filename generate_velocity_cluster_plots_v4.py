"""
OGフィルタリング済みのh5adファイルに対してvelocityを再計算し、
クラスター（OG1-4, PO1-2）で色分けしたvelocityストリームプロットを生成
"""
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

# クラスター色分けの設定
cluster_colors = {
    "OG1": "#FF0000",  # Red
    "OG2": "#FF8000",  # Orange
    "OG3": "#FFFF00",  # Yellow
    "OG4": "#00FF00",  # Green
    "PO1": "#0080FF",  # Blue
    "PO2": "#8000FF"   # Purple
}

# 処理対象のh5adファイル
h5ad_files = [
    "E17_batch_1_possorted_genome_bam_D5AOD_velocity_OG_filtered_v2.h5ad",
    "E17_batch_2_possorted_genome_bam_Z9D92_velocity_OG_filtered_v2.h5ad",
    "E17_batch_3_possorted_genome_bam_WACA8_velocity_OG_filtered_v2.h5ad"
]

for h5ad_path in h5ad_files:
    if not Path(h5ad_path).exists():
        print(f"File not found: {h5ad_path}")
        continue

    sample_name = h5ad_path.replace("_velocity_OG_filtered_v2.h5ad", "")
    print(f"\n{'='*60}")
    print(f"Processing {sample_name}")
    print('='*60)

    # Load data
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Check if cluster info exists
    if "cluster" not in adata.obs:
        print("  Warning: No cluster information found")
        continue

    # Check cluster distribution
    cluster_counts = adata.obs["cluster"].value_counts().sort_index()
    print(f"Cluster distribution:")
    for cluster, count in cluster_counts.items():
        print(f"  {cluster}: {count} cells")

    # Recompute UMAP if missing
    if "umap" not in adata.obsm:
        print("  Recomputing UMAP...")
        try:
            # Ensure we have PCA
            if "X_pca" not in adata.obsm:
                sc.pp.scale(adata, max_value=10)
                sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))

            # Compute neighbors and UMAP
            sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
            sc.tl.umap(adata)
            print("  UMAP recomputed successfully")
        except Exception as e:
            print(f"  UMAP recomputation failed: {e}")
            continue

    # Recompute velocity if missing
    if "velocity_umap" not in adata.obsm:
        print("  Recomputing velocity...")
        try:
            # Filter and normalize for velocity (scVelo 0.3.4 compatible)
            scv.pp.filter_and_normalize(adata, min_shared_counts=20)

            # Select highly variable genes if needed
            if adata.shape[1] > 2000:
                sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
                print(f"  Selected {adata.shape[1]} highly variable genes")

            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

            # Try velocity computation
            try:
                scv.tl.velocity(adata, mode="stochastic")
                print("  Stochastic velocity computed")
            except Exception as e:
                print(f"  Stochastic failed: {e}")
                try:
                    scv.tl.recover_dynamics(adata, n_jobs=4)
                    scv.tl.velocity(adata, mode="dynamical")
                    print("  Dynamical velocity computed")
                except Exception as e2:
                    print(f"  Dynamical also failed: {e2}")
                    continue

            scv.tl.velocity_graph(adata)
            print("  Velocity graph computed")

        except Exception as e:
            print(f"  Velocity computation failed: {e}")
            continue

    # Check if velocity data exists now
    has_velocity = "velocity_umap" in adata.obsm
    print(f"Has velocity_umap: {has_velocity}")

    # Create velocity stream plot with cluster colors
    try:
        fig, ax = plt.subplots(figsize=(12, 10))

        # Create custom color palette
        color_palette = [cluster_colors.get(c, "#808080") for c in sorted(adata.obs["cluster"].unique())]

        # Velocity stream plot
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="cluster",
            palette=color_palette,
            ax=ax,
            show=False,
            legend_loc="right margin",
            size=80,
            alpha=0.7,
            linewidth=1.5
        )

        ax.set_title(f"{sample_name}\nVelocity Stream (OG1-4, PO1-2)", fontsize=16, fontweight="bold")

        # Add legend with cluster colors
        handles = []
        labels = []
        for cluster in sorted(adata.obs["cluster"].unique()):
            if cluster in cluster_colors:
                handles.append(plt.Rectangle((0,0),1,1, facecolor=cluster_colors[cluster], edgecolor='black'))
                labels.append(f"{cluster} (n={cluster_counts.get(cluster, 0)})")

        ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.05, 1), title="Clusters")

        # Save plot
        out_file = OUT_DIR / f"{sample_name}_velocity_stream_clusters.png"
        fig.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out_file.name}")

        # Save updated h5ad with velocity
        adata.write(h5ad_path)
        print(f"Updated h5ad with velocity: {h5ad_path}")

    except Exception as e:
        print(f"Velocity stream plot failed: {e}")

        # Fallback: just UMAP plot with clusters
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.umap(adata, color="cluster", ax=ax, show=False, palette=color_palette, size=80)
            ax.set_title(f"{sample_name}\nUMAP Clusters (OG1-4, PO1-2)", fontsize=14, fontweight="bold")

            out_file = OUT_DIR / f"{sample_name}_umap_clusters_only.png"
            fig.savefig(out_file, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved fallback UMAP: {out_file.name}")
        except Exception as e2:
            print(f"Fallback plot also failed: {e2}")

print("\n" + "="*60)
print("All velocity plots generated.")
print("="*60)
