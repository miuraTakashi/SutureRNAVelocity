"""
OGフィルタリング済みのh5adファイルに対して、
クラスター（OG1-4, PO1-2）で色分けしたvelocityストリームプロットを生成
"""
import scanpy as sc
import scvelo as scv
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

    # Check if UMAP exists
    if "umap" not in adata.obsm:
        print("  Warning: No UMAP coordinates found")
        continue

    # Check if velocity data exists
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
