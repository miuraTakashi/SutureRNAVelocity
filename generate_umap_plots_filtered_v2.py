"""
メタデータクラスター情報を追加し、フィルタリング時にUMAPを保持。
その後UMAPプロットを生成。
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
    meta_df = pd.read_csv(f, index_col=0)

barcode_to_cluster = dict(zip(meta_df.index, meta_df['Cluster']))

h5ad_files = {
    1: "E17_batch_1_possorted_genome_bam_D5AOD_velocity.h5ad",
    2: "E17_batch_2_possorted_genome_bam_Z9D92_velocity.h5ad",
    3: "E17_batch_3_possorted_genome_bam_WACA8_velocity.h5ad"
}

def extract_and_normalize_barcode(bc_string):
    if ":" in bc_string:
        bc = bc_string.split(":")[-1]
    else:
        bc = bc_string
    bc = bc.rstrip("x")
    if bc in barcode_to_cluster:
        return bc
    if f"{bc}-1" in barcode_to_cluster:
        return f"{bc}-1"
    if f"{bc}-2" in barcode_to_cluster:
        return f"{bc}-2"
    if f"{bc}-1_1" in barcode_to_cluster:
        return f"{bc}-1_1"
    if f"{bc}-2_1" in barcode_to_cluster:
        return f"{bc}-2_1"
    return None

for batch_num, h5ad_path in h5ad_files.items():
    sample_name = h5ad_path.replace("_velocity.h5ad", "")
    print(f"\n{'='*60}")
    print(f"Processing Batch {batch_num}: {sample_name}")
    print('='*60)

    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Map clusters
    adata.obs["cluster_meta"] = [extract_and_normalize_barcode(bc) for bc in adata.obs_names]
    adata.obs["cluster"] = adata.obs["cluster_meta"].map(barcode_to_cluster)

    matched = adata.obs["cluster"].notna().sum()
    print(f"Matched cells: {matched} / {adata.shape[0]}")

    # Filter to OG/PO clusters
    clusters_to_keep = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
    mask = adata.obs["cluster"].isin(clusters_to_keep)
    
    print(f"Filtering to OG1-4, PO1-2: {mask.sum()} cells")

    # Create filtered object while preserving obsm
    adata_filtered = adata[mask].copy()
    
    # Check if UMAP exists and explicitly copy it
    if "umap" in adata.obsm:
        print(f"  UMAP shape in original: {adata.obsm['umap'].shape}")
        print(f"  UMAP shape after filter: {adata_filtered.obsm['umap'].shape}")
    else:
        print("  WARNING: No UMAP in original data")

    # Generate UMAP plot if it exists
    if "umap" in adata_filtered.obsm:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(adata_filtered, color="cluster", ax=ax, show=False, 
                   legend_loc="on data", size=100, palette="tab20")
        ax.set_title(f"{sample_name}\nOG1-4, PO1-2 (Metadata Clusters)", fontsize=14, fontweight="bold")
        
        outfile = OUT_DIR / f"{sample_name}_umap_OG_filtered.png"
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {outfile.name}")
    else:
        print("  Warning: No UMAP data to plot")

    # Distribution bar plot (always possible)
    if "cluster" in adata_filtered.obs:
        fig, ax = plt.subplots(figsize=(8, 5))
        cluster_counts = adata_filtered.obs["cluster"].value_counts().sort_index()
        cluster_counts.plot(kind="bar", ax=ax, color="steelblue")
        ax.set_title(f"{sample_name}: Cell Distribution by Cluster (OG/PO)")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Number of Cells")
        ax.tick_params(axis="x", rotation=45)
        
        outfile_dist = OUT_DIR / f"{sample_name}_cluster_distribution.png"
        fig.savefig(outfile_dist, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {outfile_dist.name}")

print("\nProcessing complete.")
