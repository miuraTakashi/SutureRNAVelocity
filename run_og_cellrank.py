"""
OG1-4 CellRank 解析
- root 固定なし（脱分化 OG4→OG2 を含む双方向の流れを許容）
- CellRank で initial state / terminal state を確率的に同定
- 各状態への吸収確率（fate probability）を計算
- PAGA トポロジーも可視化
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
import cellrank as cr

OUTDIR = "figures_og_exact"
sc.settings.figdir = OUTDIR
os.makedirs(OUTDIR, exist_ok=True)

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

OG_COLORS = {"OG1": "#e41a1c", "OG2": "#ff7f00",
             "OG3": "#4daf4a", "OG4": "#984ea3"}
OG_PALETTE = list(OG_COLORS.values())

# ============================================================
# 1. データ読み込み
# ============================================================
print("Loading E17_og_exact_velocity.h5ad ...")
adata = sc.read_h5ad("E17_og_exact_velocity.h5ad")
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4"],
)
adata.uns["og_cluster_colors"] = OG_PALETTE

print("Cluster composition:")
print(adata.obs["og_cluster"].value_counts().to_string())

# ============================================================
# 2. velocity graph の確認（root 固定なし）
# ============================================================
if "velocity_graph" not in adata.uns:
    scv.tl.velocity_graph(adata)

# velocity stream をそのまま表示（脱分化を含む流れ）
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="og_cluster",
    title="OG1-4 velocity stream (no root constraint)",
    dpi=150, save="_og_noroot_stream.png",
)

# ============================================================
# 3. PAGA で遷移確率を可視化
# ============================================================
print("\nComputing PAGA ...")
sc.tl.paga(adata, groups="og_cluster")

paga_conn = pd.DataFrame(
    adata.uns["paga"]["connectivities"].toarray(),
    index=adata.obs["og_cluster"].cat.categories,
    columns=adata.obs["og_cluster"].cat.categories,
)
print("PAGA connectivity:")
print(paga_conn.round(3).to_string())

fig, ax = plt.subplots(figsize=(5, 5))
sc.pl.paga(
    adata,
    color=["og_cluster"],
    node_size_scale=3,
    edge_width_scale=2,
    min_edge_width=1,
    threshold=0.01,
    title="PAGA: OG1-4 topology\n(edge width = transition strength)",
    ax=ax, show=False,
)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/og_paga_graph.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/og_paga_graph.png")

# ============================================================
# 4. CellRank: velocity kernel で遷移行列を構築
# ============================================================
print("\nCellRank: building VelocityKernel ...")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

# velocity (80%) + connectivity (20%) の組み合わせ
combined_kernel = 0.8 * vk + 0.2 * ck
print(f"  Kernel: {combined_kernel}")

# ============================================================
# 4.5 OG クラスター間の「分化方向」ベクトル場（粗視化）
# ============================================================
print("\nComputing coarse-grained differentiation vector field ...")

T = combined_kernel.transition_matrix  # (cells x cells) の遷移確率行列（疎行列）
clusters = adata.obs["og_cluster"].cat.categories

# 各 OG クラスターの UMAP 重心
if "X_umap" not in adata.obsm:
    raise KeyError("`adata.obsm['X_umap']` が見つかりません。UMAP 座標を含む AnnData を入力してください。")

umap = adata.obsm["X_umap"]
centroids = {}
cluster_indices = {}
for cl in clusters:
    idx = np.where(adata.obs["og_cluster"].values == cl)[0]
    cluster_indices[cl] = idx
    centroids[cl] = umap[idx].mean(axis=0)

# クラスター A→B への平均遷移確率（セル A から見た期待値）を計算
flux = pd.DataFrame(0.0, index=clusters, columns=clusters)
for a in clusters:
    idx_a = cluster_indices[a]
    if idx_a.size == 0:
        continue
    for b in clusters:
        idx_b = cluster_indices[b]
        if idx_b.size == 0:
            continue
        sub = T[idx_a, :][:, idx_b]  # A(行)→B(列)
        # 各 A 細胞について B への遷移確率を平均（A クラスター起点の期待値）
        flux.loc[a, b] = np.asarray(sub.sum(axis=1)).mean()

# 正味フラックス（A→B - B→A）に基づき、主要な分化方向のみを描画
net_flux = flux - flux.T
positive_entries = net_flux.values[net_flux.values > 0]
net_thresh = np.percentile(positive_entries, 50) if positive_entries.size > 0 else 0.0
max_net = positive_entries.max() if positive_entries.size > 0 else 1.0

fig, ax = plt.subplots(figsize=(6, 5))

# 背景に全細胞 UMAP を淡く表示
ax.scatter(umap[:, 0], umap[:, 1], s=5, c="#dddddd", alpha=0.3, linewidths=0)

# クラスター重心とラベル
for cl in clusters:
    x, y = centroids[cl]
    ax.scatter(x, y, s=80, c=OG_COLORS.get(cl, "black"), edgecolor="k", zorder=3)
    ax.text(x, y, cl, ha="center", va="center", fontsize=9, weight="bold", color="white", zorder=4)

# 正味フラックスが正でしきい値を超える向きだけ矢印として描画
for i, a in enumerate(clusters):
    for j, b in enumerate(clusters):
        if a == b:
            continue
        w = net_flux.loc[a, b]
        if w <= net_thresh:
            continue
        x0, y0 = centroids[a]
        x1, y1 = centroids[b]
        # 線幅と矢印長をフラックスに応じてスケーリング
        width = 0.5 + 3.0 * (w / max_net)
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="->",
                linewidth=width,
                color="black",
                alpha=0.8,
                shrinkA=5,
                shrinkB=5,
            ),
            zorder=2,
        )

ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title("Coarse-grained differentiation vector field\n(based on net flux between OG clusters)")
plt.tight_layout()
plt.savefig(f"{OUTDIR}/cr_diff_vector_field.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/cr_diff_vector_field.png")

# ============================================================
# 5. GPCCA で macrostate（大局的状態）を同定
# ============================================================
print("\nGPCCA: identifying macrostates ...")
estimator = cr.estimators.GPCCA(combined_kernel)

# n_states=4 (OG1〜OG4 を想定)
estimator.compute_macrostates(n_states=4, cluster_key="og_cluster")
print(f"  Macrostates identified:")
print(estimator.macrostates.value_counts().to_string())

estimator.plot_macrostates(which="all", basis="umap",
                           title="CellRank macrostates")
plt.savefig(f"{OUTDIR}/cr_macrostates.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR}/cr_macrostates.png")

# ============================================================
# 6. ループを許容したトポロジーの解釈
# ============================================================
# ここでは「必ずどこかの終末状態に吸収される」という仮定を置かず，
# 非対称（非可逆）な遷移行列に対して GPCCA で得られた macrostates の
# 構造と PAGA トポロジーを主に解釈対象とする。
# そのため terminal / initial state や absorption probability は計算しない。

# ============================================================
# 7. velocity stream を OG サブセット用ディレクトリにコピー
# ============================================================
import shutil, glob
for f in (glob.glob("figures/scvelo__og_noroot*.png")):
    dst = os.path.join(OUTDIR, os.path.basename(f).replace("scvelo__og_noroot_", "noroot_"))
    shutil.copy(f, dst)

# 保存
adata.write("E17_og_exact_velocity.h5ad")

print("\n=== DONE ===")
print(f"Figures: {OUTDIR}/")
print("Key outputs:")
print(f"  {OUTDIR}/og_paga_graph.png         -- 遷移強度グラフ")
print(f"  {OUTDIR}/cr_macrostates.png        -- CellRank 大局的状態")
print(f"  {OUTDIR}/cr_initial_terminal.png   -- 初期/終末状態")
print(f"  {OUTDIR}/cr_fate_probabilities.png -- 各状態への到達確率")
print(f"  {OUTDIR}/cr_fate_heatmap.png       -- OGクラスターごとの運命確率")
