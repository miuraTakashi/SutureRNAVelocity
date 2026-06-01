"""
単一loomから生成されたh5adファイルに対して：
1. Leidenクラスタリングを実行
2. クラスター別のUMAPプロット（ラベル付き）を生成
3. マーカー遺伝子をランキングしてセルタイプ推測に使用
4. クラスタリング結果を保存
"""
import os
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

# 処理対象のh5adファイル
h5ad_files = sorted(WORKDIR.glob("E17_batch_*_velocity.h5ad"))
if not h5ad_files:
    raise FileNotFoundError("No batch h5ad files found")

for h5ad_path in h5ad_files:
    sample_name = h5ad_path.stem.replace("_velocity", "")
    print(f"\n{'='*60}")
    print(f"Processing {sample_name}")
    print('='*60)

    # Load data
    adata = sc.read_h5ad(str(h5ad_path))
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Leiden clustering
    print("Performing Leiden clustering...")
    try:
        sc.tl.leiden(adata, resolution=1.0, key_added="leiden_cluster")
        print(f"Leiden clusters: {adata.obs['leiden_cluster'].nunique()} clusters")
        print(f"Cluster distribution:\n{adata.obs['leiden_cluster'].value_counts().sort_index()}")
    except Exception as e:
        print(f"Leiden clustering failed: {e}")
        adata.obs["leiden_cluster"] = "0"

    # Full UMAP plot with cluster labels
    if "umap" in adata.obsm:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color="leiden_cluster", ax=ax, show=False, 
                   legend_loc="on data", size=80, palette="tab20")
        ax.set_title(f"{sample_name}\nLeiden Clusters", fontsize=14, fontweight="bold")
        fig.savefig(OUT_DIR / f"{sample_name}_umap_leiden_clusters.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {sample_name}_umap_leiden_clusters.png")

    # Ranking genes and saving marker genes for cell type identification
    print("Ranking genes by cluster...")
    try:
        sc.tl.rank_genes_groups(adata, groupby="leiden_cluster", method="wilcoxon")
        marker_genes = sc.get.rank_genes_groups_df(adata, group=None)
        
        # Save marker genes
        marker_file = WORKDIR / f"{sample_name}_marker_genes.csv"
        marker_genes.to_csv(marker_file, index=False)
        print(f"Marker genes saved: {marker_file}")
        
        # Show top markers per cluster
        print("\nTop marker genes per cluster:")
        for cluster in sorted(adata.obs["leiden_cluster"].unique()):
            top_genes = marker_genes[marker_genes['group'] == cluster].head(3)['names'].tolist()
            print(f"  Cluster {cluster}: {', '.join(top_genes)}")
    except Exception as e:
        print(f"Marker gene ranking failed: {e}")

    # Save full h5ad with clusters
    out_h5ad_full = WORKDIR / f"{sample_name}_velocity_clustered.h5ad"
    adata.write(out_h5ad_full)
    print(f"Saved clustered AnnData: {out_h5ad_full}")

    # Create a summary table
    summary_file = WORKDIR / f"{sample_name}_cluster_summary.csv"
    summary_df = pd.DataFrame({
        'sample': sample_name,
        'n_clusters': adata.obs['leiden_cluster'].nunique(),
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
    }, index=[0])
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved summary: {summary_file}")

    # Velocity stream plot colored by cluster
    try:
        import scvelo as scv
        ax = scv.pl.velocity_embedding_stream(adata, basis="umap", color="leiden_cluster", 
                                               show=False, legend_loc="on data", palette="tab20")
        ax.figure.savefig(OUT_DIR / f"{sample_name}_velocity_stream_clusters.png", dpi=300, bbox_inches="tight")
        plt.close(ax.figure)
        print(f"Saved velocity stream: {sample_name}_velocity_stream_clusters.png")
    except Exception as e:
        print(f"Velocity stream plot failed: {e}")

print("\n" + "="*60)
print("All samples processed successfully.")
print("="*60)
print("\nNext steps:")
print("1. Review marker genes to identify which Leiden clusters correspond to OG/PO")
print("2. Use sc.pp.filter_genes_dispersion() or marker-based filtering for downstream analysis")
print("3. Filter cells based on cluster annotations in downstream scripts")
