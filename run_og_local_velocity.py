"""
OG1–4 + PO1–2 局所 velocity 可視化
- グリッド密度・平滑化を変えて局所の流れを確認
- 各クラスター内で velocity を分離表示
- 位相ポートレート（spliced vs unspliced）でループか単調変化かを判定
- PCA 空間での velocity（UMAP の投影歪みを除く）
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")

import scanpy as sc
import scvelo as scv

OUTDIR = "figures_og_exact"
sc.settings.figdir = OUTDIR
os.makedirs(OUTDIR, exist_ok=True)

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
adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
)
adata.uns["og_cluster_colors"] = OG_PALETTE
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# ============================================================
# 2. グリッド矢印：平滑化を変えて局所/大域の流れを比較
# ============================================================
print("\n[1] Grid arrows with different smoothing ...")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
densities = [0.5, 1.0, 2.0]  # グリッド密度（高いほど細かい）

for ax, density in zip(axes, densities):
    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="og_cluster",
        density=density,
        arrow_size=2.5,
        arrow_length=2.5,
        smooth=0.5,         # 平滑化を小さくして局所の流れを強調
        ax=ax,
        show=False,
        title=f"Grid velocity (density={density})",
    )
plt.tight_layout()
plt.savefig(f"{OUTDIR}/local_grid_velocity.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/local_grid_velocity.png")

# 平滑化パラメータを変えた比較
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
smooths = [0.3, 0.8, 1.5]

for ax, sm in zip(axes, smooths):
    scv.pl.velocity_embedding_grid(
        adata, basis="umap", color="og_cluster",
        density=1.5, smooth=sm,
        arrow_size=2, arrow_length=2,
        ax=ax, show=False,
        title=f"smooth={sm}",
    )
plt.tight_layout()
plt.savefig(f"{OUTDIR}/local_smooth_compare.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/local_smooth_compare.png")

# ============================================================
# 3. 各 OG クラスター内での局所 velocity（クラスター分離表示）
# ============================================================
print("\n[2] Per-cluster local velocity ...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()

# 全体の UMAP 座標範囲を取得
umap_all = adata.obsm["X_umap"]
x_min, x_max = umap_all[:, 0].min() - 0.5, umap_all[:, 0].max() + 0.5
y_min, y_max = umap_all[:, 1].min() - 0.5, umap_all[:, 1].max() + 0.5

for i, (og, col) in enumerate(OG_COLORS.items()):
    ax = axes[i]
    # 全細胞を背景に薄く表示
    ax.scatter(umap_all[:, 0], umap_all[:, 1],
               c="#dddddd", s=5, alpha=0.3, zorder=1)
    # 対象クラスターの細胞を強調
    mask = adata.obs["og_cluster"] == og
    umap_og = umap_all[mask]
    ax.scatter(umap_og[:, 0], umap_og[:, 1],
               c=col, s=15, alpha=0.8, label=og, zorder=2)
    # グリッド矢印（全体の velocity graph を使いつつ、対象クラスター周辺に絞る）
    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="og_cluster",
        density=2.0,
        smooth=0.4,
        arrow_size=2,
        arrow_length=2,
        ax=ax, show=False,
        title=f"Local velocity around {og}",
        alpha=0.0,   # 散布点は非表示（上で描画済み）
    )
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_title(f"Local velocity (highlighted: {og})", fontsize=11)
    ax.legend(loc="upper right", markerscale=2)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/local_per_cluster_stream.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/local_per_cluster_stream.png")

# ============================================================
# 4. 位相ポートレート（spliced vs unspliced）
#    ループ = 時計回り/反時計回り / 単調変化 = 斜め直線上の流れ
# ============================================================
print("\n[3] Phase portraits for top velocity genes ...")

# velocity confidence の高い遺伝子を選択
if "fit_likelihood" in adata.var.columns:
    top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).head(12).index.tolist()
else:
    # fit_pars から likelihood を取得
    top_genes = adata.var_names[:12].tolist()

print(f"  Top genes: {top_genes}")

# 位相ポートレート（4×3 グリッド）
fig, axes = plt.subplots(3, 4, figsize=(18, 12))
axes = axes.flatten()

for i, gene in enumerate(top_genes[:12]):
    ax = axes[i]
    if gene not in adata.var_names:
        ax.axis("off")
        continue
    try:
        scv.pl.scatter(
            adata,
            x="Ms",    # spliced（平滑化済み）
            y="Mu",    # unspliced（平滑化済み）
            color="og_cluster",
            var_names=gene,
            ax=ax, show=False,
            title=gene,
            legend_loc=None,
            alpha=0.6,
            size=20,
        )
    except Exception as e:
        ax.set_title(f"{gene}\n(error)")
        ax.axis("off")

plt.tight_layout()
plt.savefig(f"{OUTDIR}/local_phase_portraits.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/local_phase_portraits.png")

# ============================================================
# 5. PCA 空間での velocity
#    UMAP では距離が歪むため、PCA で確認するとループが消えることがある
# ============================================================
print("\n[4] Velocity in PCA space ...")

scv.pl.velocity_embedding_stream(
    adata,
    basis="pca",
    color="og_cluster",
    title="Velocity stream in PCA space\n(less distortion than UMAP)",
    dpi=150,
    save="_og_pca_stream.png",
)

scv.pl.velocity_embedding_grid(
    adata,
    basis="pca",
    color="og_cluster",
    density=1.5,
    smooth=0.5,
    arrow_size=2,
    title="Velocity grid in PCA space",
    dpi=150,
    save="_og_pca_grid.png",
)
print(f"  Velocity in PCA space saved")

# ============================================================
# 6. velocity の方向一致スコア（velocity coherence）を可視化
#    近傍細胞との velocity 方向の一致度 → ループなら低くなる
# ============================================================
print("\n[5] Velocity coherence ...")

scv.tl.velocity_confidence(adata)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
scv.pl.scatter(
    adata, basis="umap",
    color="velocity_confidence",
    color_map="RdYlGn",
    vmin=0, vmax=1,
    title="Velocity confidence\n(high=coherent, low=loop/noise)",
    ax=axes[0], show=False,
)
# クラスターごとの velocity confidence 分布
clusters = list(OG_COLORS.keys())
conf_data = [
    adata.obs[adata.obs["og_cluster"] == og]["velocity_confidence"].dropna().values
    for og in clusters
]
axes[1].boxplot(
    conf_data,
    labels=clusters,
    patch_artist=True,
    boxprops=dict(facecolor="lightblue"),
)
for patch, color in zip(
    axes[1].findobj(matplotlib.patches.PathPatch), OG_PALETTE
):
    patch.set_facecolor(color)
axes[1].set_ylabel("Velocity confidence")
axes[1].set_title(
    "Velocity confidence per cluster\n(< 0.4 suggests loop/noise)"
)
axes[1].axhline(0.4, color="red", linestyle="--", alpha=0.7, label="threshold=0.4")
axes[1].legend()

plt.tight_layout()
plt.savefig(f"{OUTDIR}/local_velocity_confidence.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/local_velocity_confidence.png")

# confidence の統計値を表示
conf_stats = adata.obs.groupby("og_cluster", observed=True)[
    "velocity_confidence"
].agg(["mean", "median", "std"]).round(3)
print("\nVelocity confidence per OG cluster:")
print(conf_stats.to_string())
print("(1.0 = perfect coherence, < 0.4 = likely loop or noise)")

# ============================================================
# コピー
# ============================================================
import shutil, glob
for f in (glob.glob("figures/scvelo__og_pca*.png")):
    dst = os.path.join(OUTDIR, os.path.basename(f).replace("scvelo__og_pca_", "pca_"))
    shutil.copy(f, dst)

print("\n=== DONE ===")
print(f"Figures saved to {OUTDIR}/:")
print("  local_grid_velocity.png      -- グリッド矢印（密度比較）")
print("  local_smooth_compare.png     -- 平滑化パラメータ比較")
print("  local_per_cluster_stream.png -- OGクラスター別ストリーム")
print("  local_phase_portraits.png    -- 位相ポートレート（ループ判定）")
print("  pca_stream.png / pca_grid.png -- PCA空間での velocity")
print("  local_velocity_confidence.png -- velocity 一貫性スコア")
