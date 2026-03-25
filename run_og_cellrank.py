"""
OG1–4 + PO1–2 CellRank 解析
- root 固定なし（双方向の流れとループを許容）
- CellRank で macrostates を同定し、PAGA トポロジーとともに解釈
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

# Batch correction 用（Harmony または他の方法）
try:
    from harmony import harmony
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False
    print("Warning: harmonyphpy not installed. Will use alternative batch correction.")

OUTDIR = "figures_og_exact"
OUTDIR_BATCH = "figures_og_exact_batch_corrected"
os.makedirs(OUTDIR_BATCH, exist_ok=True)
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
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
)
adata.uns["og_cluster_colors"] = OG_PALETTE

print("Cluster composition:")
print(adata.obs["og_cluster"].value_counts().to_string())

# ============================================================
# 1.5 Batch 情報の確認
# ============================================================
print("\n=== Batch Information ===")
if "batch" in adata.obs:
    print("Batch annotation found in adata.obs['batch']:")
    print(adata.obs["batch"].value_counts().to_string())
    batch_key = "batch"
elif "SRR" in adata.obs.columns or "run_accession" in adata.obs.columns:
    print("SRR/run_accession found, using as batch identifier.")
    if "SRR" in adata.obs.columns:
        batch_key = "SRR"
    else:
        batch_key = "run_accession"
    print(adata.obs[batch_key].value_counts().to_string())
else:
    print("No explicit batch column found. Checking metadata...")
    print(f"Available columns: {list(adata.obs.columns)}")
    print("\nNote: If no batch information, batch effect may still be present.")
    batch_key = None

# ============================================================
# 1.6 Batch Correction（Harmony または BBKNN）
# ============================================================
if batch_key:
    print(f"\n=== Applying Batch Correction (batch key: {batch_key}) ===")
    
    # 方法1: Harmony（推奨）
    if HAS_HARMONY:
        print("Using Harmony for batch correction...")
        adata_corrected = harmony(
            adata,
            key=batch_key,
            max_iter_harmony=10,
        )
        print("  Harmony batch correction completed.")
    else:
        # 方法2: BBKNN を使用（scanpy に組み込まれている）
        print("Using BBKNN for batch-aware neighbor graph...")
        import bbknn
        bbknn.bbknn(adata, batch_key=batch_key, n_pcs=50)
        adata_corrected = adata.copy()
        print("  BBKNN applied.")
    
    # Batch correction の効果を可視化
    print("\nVisualizing batch correction effect...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Before batch correction
    sc.pl.umap(
        adata,
        color=batch_key,
        ax=axes[0],
        show=False,
        title="Before Batch Correction",
    )
    
    # After batch correction (PCA + UMAP を再計算)
    sc.pp.pca(adata_corrected, n_comps=50)
    sc.pp.neighbors(adata_corrected, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata_corrected)
    
    sc.pl.umap(
        adata_corrected,
        color=batch_key,
        ax=axes[1],
        show=False,
        title="After Batch Correction",
    )
    
    plt.tight_layout()
    plt.savefig(f"{OUTDIR_BATCH}/batch_correction_effect.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {OUTDIR_BATCH}/batch_correction_effect.png")
    
    # og_cluster による可視化も保存
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.umap(
        adata,
        color="og_cluster",
        palette=OG_PALETTE,
        ax=axes[0],
        show=False,
        title="Before Batch Correction (colored by OG cluster)",
    )
    
    sc.pl.umap(
        adata_corrected,
        color="og_cluster",
        palette=OG_PALETTE,
        ax=axes[1],
        show=False,
        title="After Batch Correction (colored by OG cluster)",
    )
    
    plt.tight_layout()
    plt.savefig(f"{OUTDIR_BATCH}/og_clusters_before_after.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {OUTDIR_BATCH}/og_clusters_before_after.png")
    
    # Batch correction 後のデータを使用
    adata = adata_corrected
    print("Using batch-corrected data for downstream analysis.")
else:
    print("\nNo batch key found. Proceeding without batch correction.")
    print("Warning: If batch effect is present, results may be biased.")

# ============================================================
# 2. velocity graph の確認（root 固定なし）
# ============================================================
# Batch correction 後は velocity が無効になる可能性があるため、
# 元のデータから velocity を引き継いで、PCA/UMAP の座標のみ更新
if "velocity_graph" not in adata.uns:
    print("\nRecomputing velocity graph (batch-corrected PCA space)...")
    scv.tl.velocity_graph(adata)
else:
    print("\nVelocity graph preserved from original data.")

sc.settings.figdir = OUTDIR_BATCH

# velocity stream をそのまま表示（脱分化を含む流れ）
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="og_cluster",
    title="OG1–4 + PO1–2 velocity stream (no root constraint, batch-corrected)",
    dpi=150,
    save="_og_noroot_stream.png",
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
    title="PAGA: OG/PO topology (batch-corrected)\n(edge width = transition strength)",
    ax=ax,
    show=False,
)
plt.tight_layout()
plt.savefig(f"{OUTDIR_BATCH}/og_paga_graph.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR_BATCH}/og_paga_graph.png")

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
ax.set_title("Coarse-grained differentiation vector field (batch-corrected)\n(based on net flux between OG clusters)")
plt.tight_layout()
plt.savefig(f"{OUTDIR_BATCH}/cr_diff_vector_field.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR_BATCH}/cr_diff_vector_field.png")

# ============================================================
# 5. GPCCA で macrostate（大局的状態）を同定
#    ※可視化ラベルは常に元の OG/PO クラスター（og_cluster）に統一
# ============================================================
print("\nGPCCA: identifying macrostates ...")
estimator = cr.estimators.GPCCA(combined_kernel)

# n_states は 4–6 程度で十分。ここでは 6 とし、OG1–4 + PO1–2 を想定。
estimator.compute_macrostates(n_states=6, cluster_key="og_cluster")
print("  Macrostates identified (raw labels):")
print(estimator.macrostates.value_counts().to_string())

# 各 macrostate を、その中で最多の og_cluster にマッピング
macro = estimator.macrostates
major_map = {}
for state in macro.cat.categories:
    idx = macro[macro == state].index
    if len(idx) == 0:
        continue
    majority = (
        adata.obs.loc[idx, "og_cluster"]
        .value_counts()
        .idxmax()
    )
    major_map[state] = majority

adata.obs["macro_major"] = macro.replace(major_map)
adata.obs["macro_major"] = pd.Categorical(
    adata.obs["macro_major"],
    categories=adata.obs["og_cluster"].cat.categories,
)

print("\n  Macrostate → og_cluster mapping used for plots:")
for k, v in major_map.items():
    print(f"    {k} -> {v}")

# 可視化は常に OG/PO ラベルで統一
fig, ax = plt.subplots(figsize=(5, 5))
sc.pl.umap(
    adata,
    color="macro_major",
    palette=OG_PALETTE,
    ax=ax,
    show=False,
    title="CellRank macrostates (batch-corrected, colored by OG/PO cluster)",
)
plt.tight_layout()
plt.savefig(f"{OUTDIR_BATCH}/cr_macrostates.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {OUTDIR_BATCH}/cr_macrostates.png (labels unified to OG/PO)")

# ============================================================
# 6. ループを許容したトポロジーの解釈
# ============================================================
# ここでは「必ずどこかの終末状態に吸収される」という仮定を置かず，
# 非対称（非可逆）な遷移行列に対して GPCCA で得られた macrostates の
# 構造と PAGA トポロジーを主に解釈対象とする。
# そのため terminal / initial state や absorption probability は計算しない。

# ============================================================
# 6.5 Batch correction の統計サマリー
# ============================================================
print("\n=== Summary Statistics ===")
print(f"\nOG/PO cluster distribution (after batch correction):")
print(adata.obs["og_cluster"].value_counts().to_string())

# PAGA 接続性の比較（batch 補正前後）
print(f"\nPAGA connectivity (after batch correction):")
print(paga_conn.round(3).to_string())

# ============================================================
# 7. 結果の保存
# ============================================================
print("\nSaving results...")
adata.write(f"{OUTDIR_BATCH}/E17_og_exact_velocity_batch_corrected.h5ad")

print(f"\n=== BATCH CORRECTION ANALYSIS COMPLETE ===")
print(f"\nBatch-corrected results saved to: {OUTDIR_BATCH}/")
print("\nKey outputs:")
print(f"  {OUTDIR_BATCH}/batch_correction_effect.png")
print(f"  {OUTDIR_BATCH}/og_clusters_before_after.png")
print(f"  {OUTDIR_BATCH}/og_paga_graph.png         -- 遷移強度グラフ（バッチ補正済み）")
print(f"  {OUTDIR_BATCH}/cr_macrostates.png        -- CellRank 大局的状態（バッチ補正済み）")
print(f"  {OUTDIR_BATCH}/cr_diff_vector_field.png  -- 分化方向ベクトル場")
print(f"\nComparison with original (non-corrected):")
print(f"  {OUTDIR}/ (original directory)")
print(f"\nBatch-corrected data saved to:")
print(f"  {OUTDIR_BATCH}/E17_og_exact_velocity_batch_corrected.h5ad")
