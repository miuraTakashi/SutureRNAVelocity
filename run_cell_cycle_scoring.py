"""
OG/PO クラスターの cell cycle scoring を計算し、可視化するスクリプト

- E17_og_exact_velocity.h5ad を読み込み
- OG1-4, PO1, PO2 のみを対象とする
- Scanpy の cell cycle scoring を実行
- クラスターごとの score をボックス/バー/UMAP で可視化
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

OUTDIR = "figures_cell_cycle"
os.makedirs(OUTDIR, exist_ok=True)
adata_path = "E17_og_exact_velocity.h5ad"
if not os.path.exists(adata_path):
    raise FileNotFoundError(f"{adata_path} not found. Please run preprocessing first.")

print(f"Loading {adata_path}...")
adata = sc.read_h5ad(adata_path)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# 対象クラスターを絞る
target_clusters = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=target_clusters,
)
adata = adata[adata.obs["og_cluster"].isin(target_clusters)].copy()
print(f"Filtered to {adata.n_obs} cells in {len(target_clusters)} clusters.")

# UMAP がない場合は再計算
if "X_umap" not in adata.obsm:
    print("Computing UMAP...")
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)

# Cell cycle gene lists for score_genes_cell_cycle
s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN",
    "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
    "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8"
]
g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
    "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
    "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA3", "CDC45", "CDC6", "CDCA2",
    "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE",
    "CENPA"
]

# Gene names がすべて大文字でなければ変換
adata.var_names_make_unique()
adata.var_names = [g.upper() for g in adata.var_names]

print("Running cell cycle scoring...
Checking gene list against dataset...")
adata.var_names_make_unique()
adata.var_names = [g.upper() for g in adata.var_names]

valid_s = [g for g in s_genes if g in adata.var_names]
valid_g2m = [g for g in g2m_genes if g in adata.var_names]
print(f"Valid S-phase genes: {len(valid_s)} / {len(s_genes)}")
print(f"Valid G2M-phase genes: {len(valid_g2m)} / {len(g2m_genes)}")
if len(valid_s) == 0 and len(valid_g2m) == 0:
    msg = (
        "No canonical cell cycle genes were found in the dataset. "
        "Cell cycle scoring cannot be computed on this gene set."
    )
    print(msg)
    with open(os.path.join(OUTDIR, "cell_cycle_error.txt"), "w") as f:
        f.write(msg + "\n")
    raise SystemExit(msg)

sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=valid_s,
    g2m_genes=valid_g2m,
)

# 各クラスターの平均値を計算
cluster_scores = adata.obs.groupby("og_cluster")[ ["S_score", "G2M_score"] ].mean()
cluster_scores["cell_cycle_score"] = cluster_scores[["S_score", "G2M_score"]].max(axis=1)
print(cluster_scores.round(4))
cluster_scores.to_csv(os.path.join(OUTDIR, "cluster_cell_cycle_scores.csv"))

# UMAP プロット
sc.settings.figdir = OUTDIR
sc.pl.umap(
    adata,
    color=["S_score", "G2M_score", "phase"],
    cmap="viridis",
    save="_cell_cycle_scores.png",
    show=False,
)
# scanpy は自動的に OUTDIR に画像を保存するため、追加のコピーは不要

# クラスターごとのボックスプロット
plot_data = adata.obs.reset_index()
fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
fig.suptitle("Cell cycle scores by OG/PO cluster")
sns.boxplot(
    data=plot_data,
    x="og_cluster",
    y="S_score",
    order=target_clusters,
    palette="Set2",
    ax=axes[0],
)
axes[0].set_title("S phase score")
axes[0].set_xlabel("")
axes[0].set_ylabel("S_score")
sns.boxplot(
    data=plot_data,
    x="og_cluster",
    y="G2M_score",
    order=target_clusters,
    palette="Set2",
    ax=axes[1],
)
axes[1].set_title("G2M phase score")
axes[1].set_xlabel("")
axes[1].set_ylabel("G2M_score")
plt.savefig(os.path.join(OUTDIR, "cell_cycle_boxplots.png"), dpi=200)
plt.close()

# クラスターごとの平均値棒グラフ
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
cluster_scores.plot(kind="bar", ax=ax, color=["#4c72b0", "#dd8452", "#55a868"])
ax.set_ylabel("Mean score")
ax.set_title("Mean cell cycle scores per cluster")
ax.legend(title="Score")
plt.xticks(rotation=0)
plt.savefig(os.path.join(OUTDIR, "cluster_mean_scores.png"), dpi=200)
plt.close()

print(f"Saved cell cycle visualizations to {OUTDIR}")
print("Done.")
