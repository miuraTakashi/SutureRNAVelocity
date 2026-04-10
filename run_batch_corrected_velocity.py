"""
バッチ補正を行った上でvelocityグラフを計算するスクリプト

- E17_og_exact_velocity.h5ad を読み込み
- バッチ補正（HarmonyまたはBBKNN）を適用
- PCA/UMAPを再計算
- Velocityグラフを計算
- ストリームプロットを保存
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

import scanpy as sc
import scvelo as scv

# Batch correction 用
try:
    from harmony import harmonize
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False
    print("Warning: harmonypy not installed. Will use BBKNN.")

OUTDIR = "figures_batch_corrected_velocity"
os.makedirs(OUTDIR, exist_ok=True)
sc.settings.figdir = OUTDIR

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

OG_COLORS = {
    "OG1": "#e41a1c",
    "OG2": "#ff7f00",
    "OG3": "#4daf4a",
    "OG4": "#984ea3",
    "PO1": "#377eb8",
    "PO2": "#a65628",
}
OG_PALETTE = list(OG_COLORS.values())

# ============================================================
# 1. データ読み込み
# ============================================================
print("Loading E17_og_exact_velocity.h5ad ...")
adata = sc.read_h5ad("E17_og_exact_velocity.h5ad")
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
)
adata.uns["og_cluster_colors"] = OG_PALETTE

print("Cluster composition:")
print(adata.obs["og_cluster"].value_counts().to_string())

# ============================================================
# 2. Batch 情報の確認
# ============================================================
print("\n=== Batch Information ===")
batch_key = None
if "batch" in adata.obs:
    batch_key = "batch"
elif "SRR" in adata.obs.columns:
    batch_key = "SRR"
elif "run_accession" in adata.obs.columns:
    batch_key = "run_accession"
else:
    print("No batch column found. Proceeding without batch correction.")
    print("Note: If batch effect exists, add batch annotation to adata.obs['batch'].")

if batch_key:
    print(f"Using {batch_key} as batch key.")
    print(adata.obs[batch_key].value_counts().to_string())

# ============================================================
# 3. Batch Correction
# ============================================================
if batch_key:
    print("\n=== Applying Batch Correction ===")

    if HAS_HARMONY:
        print("Using Harmony...")
        adata_corrected = harmonize(adata, key=batch_key, max_iter_harmony=10)
    else:
        print("Using BBKNN...")
        import bbknn
        bbknn.bbknn(adata, batch_key=batch_key, n_pcs=50)
        adata_corrected = adata.copy()

    # PCA/UMAP再計算
    sc.pp.pca(adata_corrected, n_comps=50)
    sc.pp.neighbors(adata_corrected, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata_corrected)

    # 可視化
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(adata, color=batch_key, ax=axes[0], show=False, title="Before Batch Correction")
    sc.pl.umap(adata_corrected, color=batch_key, ax=axes[1], show=False, title="After Batch Correction")
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/batch_correction_effect.png", dpi=150)
    plt.close()

    adata = adata_corrected
    print("Batch correction applied.")
else:
    print("No batch correction applied.")

# ============================================================
# 4. Velocity Graph 計算
# ============================================================
print("\n=== Computing Velocity Graph ===")
scv.tl.velocity_graph(adata)

# ============================================================
# 5. 可視化
# ============================================================
print("Generating velocity stream plot...")
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="og_cluster",
    palette=OG_PALETTE,
    title="Velocity Stream (Batch Corrected)",
    save="_velocity_stream_batch_corrected.png",
)

print(f"Results saved to {OUTDIR}")
print("Done.")