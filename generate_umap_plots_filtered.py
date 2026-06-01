"""
メタデータクラスター付きのUMAPプロットを生成（Velocityストリームなしの簡単版）
"""
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

h5ad_files = [
    "E17_batch_1_possorted_genome_bam_D5AOD_velocity_OG_filtered.h5ad",
    "E17_batch_2_possorted_genome_bam_Z9D92_velocity_OG_filtered.h5ad",
    "E17_batch_3_possorted_genome_bam_WACA8_velocity_OG_filtered.h5ad"
]

for h5ad_path in h5ad_files:
    if not Path(h5ad_path).exists():
        print(f"File not found: {h5ad_path}")
        continue

    sample_name = h5ad_path.replace("_velocity_OG_filtered.h5ad", "")
    print(f"Processing {sample_name}...")

    # Load filtered h5ad
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Shape: {adata.shape}")

    # UMAP plot with cluster labels
    if "umap" in adata.obsm and "cluster" in adata.obs:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color="cluster", ax=ax, show=False, 
                   legend_loc="on data", size=100, palette="tab20")
        ax.set_title(f"{sample_name}\nOG1-4, PO1-2 (Metadata Clusters)", fontsize=14, fontweight="bold")
        
        outfile = OUT_DIR / f"{sample_name}_umap_OG_filtered.png"
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outfile.name}")

        # Create cluster distribution plot
        fig, ax = plt.subplots(figsize=(8, 5))
        cluster_counts = adata.obs["cluster"].value_counts().sort_index()
        cluster_counts.plot(kind="bar", ax=ax, color="steelblue")
        ax.set_title(f"{sample_name}: Cell Distribution by Cluster")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Number of Cells")
        ax.tick_params(axis="x", rotation=45)
        
        outfile_dist = OUT_DIR / f"{sample_name}_cluster_distribution.png"
        fig.savefig(outfile_dist, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outfile_dist.name}")
    else:
        print(f"  Warning: UMAP or cluster data missing")

print("\nAll plots generated successfully.")
