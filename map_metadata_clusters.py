"""
GSE163693メタデータを参照して、Loomから生成されたh5adにクラスター情報を追加。
その後、OG1-4, PO1, PO2でフィルタリングしてプロットを生成。
"""
import gzip
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

# Load metadata
print("Loading GEO metadata...")
metadata_path = "GSE163693_E15_17_composite_metadata.csv.gz"
with gzip.open(metadata_path, "rt") as f:
    meta_df = pd.read_csv(f, index_col=0)  # Use first column as index (barcode)

print(f"Metadata shape: {meta_df.shape}")
print(f"Cluster categories: {meta_df['Cluster'].unique()}")

# Create barcode -> cluster mapping
barcode_to_cluster = dict(zip(meta_df.index, meta_df['Cluster']))
print(f"\nBarcode mapping dictionary created: {len(barcode_to_cluster)} entries")

# Process individual h5ad files
h5ad_files = {
    1: "E17_batch_1_possorted_genome_bam_D5AOD_velocity.h5ad",
    2: "E17_batch_2_possorted_genome_bam_Z9D92_velocity.h5ad",
    3: "E17_batch_3_possorted_genome_bam_WACA8_velocity.h5ad"
}

for batch_num, h5ad_path in h5ad_files.items():
    sample_name = h5ad_path.replace("_velocity.h5ad", "")
    print(f"\n{'='*60}")
    print(f"Processing Batch {batch_num}: {sample_name}")
    print('='*60)

    # Load h5ad
    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Extract barcode and map to cluster
    # Format: E17_batch_1_...:AAAACTGGx -> extract AAACCTGG and add -1 suffix
    def extract_and_normalize_barcode(bc_string):
        if ":" in bc_string:
            bc = bc_string.split(":")[-1]  # Get part after ":"
        else:
            bc = bc_string
        bc = bc.rstrip("x")  # Remove trailing x
        # Add standard 10x suffix: -1 (default) or -2
        # Check metadata for exact match first
        if bc in barcode_to_cluster:
            return bc
        # Try with -1 suffix
        if f"{bc}-1" in barcode_to_cluster:
            return f"{bc}-1"
        # Try with -2 suffix
        if f"{bc}-2" in barcode_to_cluster:
            return f"{bc}-2"
        # Try with _1 suffix (from metadata)
        if f"{bc}-1_1" in barcode_to_cluster:
            return f"{bc}-1_1"
        if f"{bc}-2_1" in barcode_to_cluster:
            return f"{bc}-2_1"
        return None

    # Map barcodes to clusters
    adata.obs["cluster_meta"] = [extract_and_normalize_barcode(bc) for bc in adata.obs_names]
    adata.obs["cluster"] = adata.obs["cluster_meta"].map(barcode_to_cluster)

    # Check mapping success
    matched = adata.obs["cluster"].notna().sum()
    print(f"Matched cells: {matched} / {adata.shape[0]} ({100*matched/adata.shape[0]:.1f}%)")

    if matched > 0:
        print(f"\nCluster distribution (matched cells):")
        matched_clusters = adata.obs[adata.obs["cluster"].notna()]["cluster"]
        print(matched_clusters.value_counts())

    # Full UMAP plot with cluster labels (if matched)
    if "umap" in adata.obsm and matched > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata, color="cluster", ax=ax, show=False, 
                   legend_loc="on data", size=80, palette="tab20")
        ax.set_title(f"{sample_name}\nClusters (from metadata)", fontsize=14, fontweight="bold")
        fig.savefig(OUT_DIR / f"{sample_name}_umap_metadata_cluster.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {sample_name}_umap_metadata_cluster.png")

    # Filter to OG1-4, PO1, PO2
    clusters_to_keep = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
    mask = adata.obs["cluster"].isin(clusters_to_keep)
    adata_filtered = adata[mask].copy()

    print(f"\nAfter filtering to {clusters_to_keep}:")
    print(f"  Cells: {adata_filtered.shape[0]} (from {adata.shape[0]})")
    if adata_filtered.shape[0] > 0:
        print(f"  Distribution:\n{adata_filtered.obs['cluster'].value_counts()}")

    # Save filtered h5ad
    out_h5ad_filtered = WORKDIR / f"{sample_name}_velocity_OG_filtered.h5ad"
    adata_filtered.write(out_h5ad_filtered)
    print(f"Saved filtered AnnData: {out_h5ad_filtered}")

    # Filtered UMAP plot
    if "umap" in adata_filtered.obsm and adata_filtered.shape[0] > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata_filtered, color="cluster", ax=ax, show=False, 
                   legend_loc="on data", size=100, palette="tab20")
        ax.set_title(f"{sample_name} (OG1-4, PO1-2)\nClusters", fontsize=14, fontweight="bold")
        fig.savefig(OUT_DIR / f"{sample_name}_umap_OG_filtered.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {sample_name}_umap_OG_filtered.png")

    # Velocity stream plot for filtered data
    if adata_filtered.shape[0] > 0:
        try:
            import scvelo as scv
            ax = scv.pl.velocity_embedding_stream(adata_filtered, basis="umap", color="cluster", 
                                                   show=False, legend_loc="on data", palette="tab20")
            ax.figure.savefig(OUT_DIR / f"{sample_name}_velocity_stream_OG_filtered.png", dpi=300, bbox_inches="tight")
            plt.close(ax.figure)
            print(f"Saved: {sample_name}_velocity_stream_OG_filtered.png")
        except Exception as e:
            print(f"Velocity stream plot failed: {e}")

print("\n" + "="*60)
print("All samples processed successfully.")
print("="*60)
