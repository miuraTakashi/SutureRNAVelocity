"""
バッチ別のUMAPプロットを生成
"""

import scanpy as sc
import matplotlib.pyplot as plt
import os

# 出力ディレクトリ
outdir = "figures_batch_umap"
os.makedirs(outdir, exist_ok=True)

# データ読み込み
adata = sc.read_h5ad('E17_og_integrated_harmony_filtered.h5ad')

print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
print(f"Batch categories: {adata.obs['batch'].unique()}")

# UMAPプロット: batch別
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(adata, color='batch', ax=ax, show=False, legend_loc='on data')
plt.savefig(f"{outdir}/umap_by_batch.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {outdir}/umap_by_batch.png")

# 比較のため、og_cluster別も生成
if 'og_cluster' in adata.obs:
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='og_cluster', ax=ax, show=False, legend_loc='on data')
    plt.savefig(f"{outdir}/umap_by_og_cluster.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outdir}/umap_by_og_cluster.png")

print("Batch UMAP plots generated.")