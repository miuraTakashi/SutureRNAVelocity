"""
クラスタ塗り分けUMAP上にvelocity stream plotを重ね、クラスター名の凡例を表示する。
"""
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

cluster_colors = {
    "OG1": "#FF0000",
    "OG2": "#FF8000",
    "OG3": "#FFFF00",
    "OG4": "#00FF00",
    "PO1": "#0080FF",
    "PO2": "#8000FF"
}

h5ad_files = [
    "E17_batch_1_possorted_genome_bam_D5AOD_velocity_OG_filtered_v2.h5ad",
    "E17_batch_2_possorted_genome_bam_Z9D92_velocity_OG_filtered_v2.h5ad",
    "E17_batch_3_possorted_genome_bam_WACA8_velocity_OG_filtered_v2.h5ad"
]

for h5ad_path in h5ad_files:
    if not Path(h5ad_path).exists():
        print(f"File not found: {h5ad_path}")
        continue

    sample_name = Path(h5ad_path).stem.replace("_velocity_OG_filtered_v2", "")
    print(f"\n{'='*60}")
    print(f"Processing {sample_name}")
    print('='*60)

    adata = sc.read_h5ad(h5ad_path)
    if "cluster" not in adata.obs:
        print("  No cluster info found. Skipping.")
        continue

    if "X_umap" not in adata.obsm:
        print("  Recomputing UMAP...")
        if "X_pca" not in adata.obsm:
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))
        sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
        sc.tl.umap(adata)
        print("  UMAP recomputed.")

    velocity_graph_bad = False
    if "velocity_graph" in adata.uns and adata.uns["velocity_graph"].shape != (adata.n_obs, adata.n_obs):
        velocity_graph_bad = True
        print(f"  Existing velocity_graph has wrong shape {adata.uns['velocity_graph'].shape} for {adata.n_obs} cells.")

    if "velocity_graph" not in adata.uns or adata.uns.get("velocity_graph") is None or velocity_graph_bad:
        print("  Recomputing velocity graph...")
        try:
            for key in ["velocity_graph", "velocity_graph_neg", "velocity_params"]:
                adata.uns.pop(key, None)
            adata.obsm.pop("velocity_umap", None)
            if "X_pca" not in adata.obsm:
                sc.pp.scale(adata, max_value=10)
                sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))
            sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
            scv.tl.velocity_graph(adata)
            scv.tl.velocity_embedding(adata, basis="umap")
            print("  Velocity graph and embedding recomputed.")
        except Exception as e:
            print(f"  Velocity graph failed: {e}")
            print("  Skipping stream overlay for this sample.")
            continue

    unique_clusters = sorted(adata.obs["cluster"].unique(), key=lambda x: str(x))
    palette = [cluster_colors.get(cluster, "#999999") for cluster in unique_clusters]

    fig, ax = plt.subplots(figsize=(12, 10))

    # Cluster background scatter
    for cluster, color in zip(unique_clusters, palette):
        mask = adata.obs["cluster"] == cluster
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            s=30,
            c=color,
            alpha=0.45,
            edgecolors="none",
            label=cluster
        )

    # Overlay velocity stream
    try:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="cluster",
            palette=palette,
            ax=ax,
            show=False,
            legend_loc="none",
            linewidth=1.5,
            alpha=0.9,
            size=50
        )
    except Exception as e:
        print(f"  Stream plot overlay failed: {e}")
        print("  Saving cluster-only plot instead.")

    ax.set_title(f"{sample_name} - Cluster + Velocity Stream", fontsize=16, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.grid(False)

    # Add legend for clusters
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cluster_colors.get(cluster, '#999999'), markersize=10) for cluster in unique_clusters]
    ax.legend(handles, unique_clusters, title="Clusters", bbox_to_anchor=(1.05, 1), loc="upper left")

    out_file = OUT_DIR / f"{sample_name}_cluster_stream_overlay.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved overlay plot: {out_file.name}")

print("\nAll cluster stream overlay plots generated.")
