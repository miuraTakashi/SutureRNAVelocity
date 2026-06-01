"""
メタデータクラスター情報を追加し、OG1-4/PO1-2でフィルタリング。
クラスター分布図とサマリーを生成。
"""
import gzip
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

# Load metadata
print("Loading metadata...")
metadata_path = "GSE163693_E15_17_composite_metadata.csv.gz"
with gzip.open(metadata_path, "rt") as f:
    meta_df = pd.read_csv(f, index_col=0)

barcode_to_cluster = dict(zip(meta_df.index, meta_df['Cluster']))

def extract_barcode(bc_string):
    if ":" in bc_string:
        bc = bc_string.split(":")[-1]
    else:
        bc = bc_string
    bc = bc.rstrip("x")
    if bc in barcode_to_cluster:
        return bc
    for suffix in ["-1", "-2", "-1_1", "-2_1"]:
        if f"{bc}{suffix}" in barcode_to_cluster:
            return f"{bc}{suffix}"
    return None

# Process each loom h5ad
h5ad_files = {
    1: "E17_batch_1_possorted_genome_bam_D5AOD_velocity.h5ad",
    2: "E17_batch_2_possorted_genome_bam_Z9D92_velocity.h5ad",
    3: "E17_batch_3_possorted_genome_bam_WACA8_velocity.h5ad"
}

summary_records = []

for batch_num, h5ad_path in h5ad_files.items():
    sample_name = h5ad_path.replace("_velocity.h5ad", "")
    print(f"\nProcessing {sample_name}...")

    adata = sc.read_h5ad(h5ad_path)
    
    # Map clusters
    adata.obs["barcode_extracted"] = [extract_barcode(bc) for bc in adata.obs_names]
    adata.obs["cluster"] = adata.obs["barcode_extracted"].map(barcode_to_cluster)

    matched = adata.obs["cluster"].notna().sum()
    print(f"  Total cells: {adata.shape[0]}")
    print(f"  Matched with metadata: {matched} ({100*matched/adata.shape[0]:.1f}%)")

    # Filter to OG/PO
    clusters_to_keep = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
    mask = adata.obs["cluster"].isin(clusters_to_keep)
    og_po_cells = mask.sum()
    
    print(f"  OG/PO cells: {og_po_cells}")
    
    if og_po_cells > 0:
        adata_og = adata[mask].copy()
        
        # Save summary
        cluster_dist = adata_og.obs["cluster"].value_counts().to_dict()
        summary_records.append({
            "Sample": sample_name,
            "Total_Cells": adata.shape[0],
            "Matched_Cells": matched,
            "OG_PO_Cells": og_po_cells,
            "OG1": cluster_dist.get("OG1", 0),
            "OG2": cluster_dist.get("OG2", 0),
            "OG3": cluster_dist.get("OG3", 0),
            "OG4": cluster_dist.get("OG4", 0),
            "PO1": cluster_dist.get("PO1", 0),
            "PO2": cluster_dist.get("PO2", 0),
        })
        
        # Create distribution plot
        fig, ax = plt.subplots(figsize=(8, 5))
        cluster_counts = adata_og.obs["cluster"].value_counts().sort_index()
        colors = {"OG1": "red", "OG2": "orange", "OG3": "gold", "OG4": "yellow", "PO1": "green", "PO2": "blue"}
        bar_colors = [colors.get(c, "gray") for c in cluster_counts.index]
        cluster_counts.plot(kind="bar", ax=ax, color=bar_colors)
        ax.set_title(f"{sample_name}: Cell Distribution (OG1-4, PO1-2)", fontsize=12, fontweight="bold")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Number of Cells")
        ax.tick_params(axis="x", rotation=45)
        
        out_file = OUT_DIR / f"{sample_name}_cluster_distribution_OG.png"
        fig.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out_file.name}")
        
        # Save filtered h5ad with cluster info
        out_h5ad = WORKDIR / f"{sample_name}_velocity_OG_filtered_v2.h5ad"
        adata_og.write(out_h5ad)
        print(f"  Saved filtered h5ad: {out_h5ad.name}")

# Create summary table
if summary_records:
    summary_df = pd.DataFrame(summary_records)
    summary_file = WORKDIR / "loom_single_analysis_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSummary table saved: {summary_file}")
    print(summary_df.to_string(index=False))

print("\nComplete.")
