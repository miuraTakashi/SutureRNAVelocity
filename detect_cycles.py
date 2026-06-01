"""
サイクル検出スクリプト for SutureRNAVelocity CellRank 解析
遷移確率行列からサイクルを検出し、OG/PO クラスター単位で分析
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

warnings.filterwarnings("ignore")

import scanpy as sc
import scvelo as scv
import cellrank as cr

# Batch correction 用
try:
    from harmony import harmony
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

OUTDIR = "figures_og_exact_cycle_analysis"
os.makedirs(OUTDIR, exist_ok=True)

OG_COLORS = {
    "OG1": "#e41a1c",
    "OG2": "#ff7f00",
    "OG3": "#4daf4a",
    "OG4": "#984ea3",
    "PO1": "#377eb8",
    "PO2": "#a65628",
}

# ============================================================
# 1. データ読み込みと CellRank 解析の再現
# ============================================================
print("Loading data and reproducing CellRank analysis...")

# バッチ補正版のデータがあるか確認
batch_corrected_file = "figures_og_exact_batch_corrected/E17_og_exact_velocity_batch_corrected.h5ad"
original_file = "E17_og_exact_velocity.h5ad"

if os.path.exists(batch_corrected_file):
    print(f"Using batch-corrected data: {batch_corrected_file}")
    adata = sc.read_h5ad(batch_corrected_file)
    analysis_type = "batch_corrected"
else:
    print(f"Using original data: {original_file}")
    adata = sc.read_h5ad(original_file)
    analysis_type = "original"

print(f"Data: {adata.n_obs} cells x {adata.n_vars} genes")

# OG/PO クラスター設定
adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
)

# ============================================================
# 2. CellRank カーネルの再構築
# ============================================================
print("\nRebuilding CellRank kernels...")

# Velocity kernel
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

# Connectivity kernel
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()

# 統合カーネル (80% velocity + 20% connectivity)
combined_kernel = 0.8 * vk + 0.2 * ck
print(f"Combined kernel: {combined_kernel}")

# ============================================================
# 3. クラスター間遷移確率行列の作成
# ============================================================
print("\nComputing cluster-level transition matrix...")

clusters = adata.obs["og_cluster"].cat.categories
n_clusters = len(clusters)

# クラスター間遷移確率行列 (cluster x cluster)
cluster_transition = pd.DataFrame(0.0, index=clusters, columns=clusters)

# 各クラスターの細胞インデックスを取得
cluster_indices = {}
for cl in clusters:
    cluster_indices[cl] = np.where(adata.obs["og_cluster"] == cl)[0]

# 細胞間遷移行列を取得
T = combined_kernel.transition_matrix

# クラスター間遷移を計算
for i, source_cl in enumerate(clusters):
    source_idx = cluster_indices[source_cl]
    if len(source_idx) == 0:
        continue

    for j, target_cl in enumerate(clusters):
        target_idx = cluster_indices[target_cl]
        if len(target_idx) == 0:
            continue

        # sourceクラスタからtargetクラスタへの遷移確率の平均
        sub_matrix = T[source_idx, :][:, target_idx]
        if sub_matrix.size > 0:
            cluster_transition.loc[source_cl, target_cl] = np.asarray(sub_matrix.mean())

print("Cluster transition matrix:")
print(cluster_transition.round(3))

# ============================================================
# 4. NetworkX グラフの作成とサイクル検出
# ============================================================
print("\nCreating NetworkX graph and detecting cycles...")

# NetworkX グラフ作成
G = nx.DiGraph()  # 有向グラフ

# ノード追加
for cl in clusters:
    G.add_node(cl, color=OG_COLORS[cl])

# エッジ追加 (遷移確率が閾値以上の場合)
transition_threshold = 0.01  # 1%以上の遷移確率のみエッジとする
for source in clusters:
    for target in clusters:
        prob = cluster_transition.loc[source, target]
        if prob > transition_threshold:
            G.add_edge(source, target, weight=prob, probability=prob)

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ============================================================
# 5. サイクル検出
# ============================================================
print("\nDetecting cycles...")

# 全てのサイクルを検出
cycles = list(nx.simple_cycles(G))
print(f"Found {len(cycles)} cycles")

# サイクル情報を整理
cycle_info = []
for i, cycle in enumerate(cycles):
    cycle_length = len(cycle)

    # サイクル内の遷移確率の積（循環確率）
    cycle_prob = 1.0
    for j in range(cycle_length):
        source = cycle[j]
        target = cycle[(j + 1) % cycle_length]
        prob = cluster_transition.loc[source, target]
        cycle_prob *= prob

    cycle_info.append({
        'cycle_id': i + 1,
        'cycle': cycle,
        'length': cycle_length,
        'cycle_probability': cycle_prob,
        'states': ', '.join(cycle)
    })

# DataFrameに変換
cycle_df = pd.DataFrame(cycle_info)
print("\nDetected cycles:")
if len(cycle_df) > 0:
    print(cycle_df.to_string(index=False))
else:
    print("No cycles detected.")

# ============================================================
# 6. 強連結成分分析 (Strongly Connected Components)
# ============================================================
print("\nAnalyzing strongly connected components...")

# 強連結成分を検出
scc = list(nx.strongly_connected_components(G))
print(f"Found {len(scc)} strongly connected components")

scc_info = []
for i, component in enumerate(scc):
    component_size = len(component)
    component_list = sorted(list(component))

    # 成分内の遷移確率の平均
    internal_probs = []
    for source in component:
        for target in component:
            if source != target:
                prob = cluster_transition.loc[source, target]
                internal_probs.append(prob)

    avg_internal_prob = np.mean(internal_probs) if internal_probs else 0.0

    scc_info.append({
        'scc_id': i + 1,
        'size': component_size,
        'states': ', '.join(component_list),
        'avg_internal_probability': avg_internal_prob,
        'is_cycle': component_size > 1  # サイズ>1はサイクル成分の可能性
    })

scc_df = pd.DataFrame(scc_info)
print("\nStrongly connected components:")
print(scc_df.to_string(index=False))

# ============================================================
# 7. グラフ可視化
# ============================================================
print("\nCreating visualization...")

plt.figure(figsize=(12, 8))

# ノード位置（spring layout）
pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

# ノード描画
node_colors = [OG_COLORS[node] for node in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=2000, alpha=0.8)

# エッジ描画（太さ = 遷移確率）
edges = G.edges()
edge_weights = [G[u][v]['weight'] * 10 for u, v in edges]  # スケーリング
nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.6, edge_color='gray',
                      arrows=True, arrowsize=20, connectionstyle='arc3,rad=0.1')

# ノードラベル
nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold', font_color='white')

# サイクルをハイライト（赤色）
if cycles:
    cycle_edges = []
    for cycle in cycles:
        for j in range(len(cycle)):
            source = cycle[j]
            target = cycle[(j + 1) % len(cycle)]
            if G.has_edge(source, target):
                cycle_edges.append((source, target))

    nx.draw_networkx_edges(G, pos, edgelist=cycle_edges, width=3,
                          edge_color='red', alpha=0.8, arrows=True,
                          arrowsize=25, connectionstyle='arc3,rad=0.1')

plt.title(f'Cluster Transition Graph with Cycles\n({analysis_type} data)', fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/cluster_transition_cycles_{analysis_type}.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {OUTDIR}/cluster_transition_cycles_{analysis_type}.png")

# ============================================================
# 8. 結果の保存
# ============================================================
print("\nSaving results...")

# サイクル情報をCSVで保存
if len(cycle_df) > 0:
    cycle_df.to_csv(f"{OUTDIR}/detected_cycles_{analysis_type}.csv", index=False)
    print(f"Saved: {OUTDIR}/detected_cycles_{analysis_type}.csv")

# SCC情報をCSVで保存
scc_df.to_csv(f"{OUTDIR}/strongly_connected_components_{analysis_type}.csv", index=False)
print(f"Saved: {OUTDIR}/strongly_connected_components_{analysis_type}.csv")

# クラスター間遷移行列を保存
cluster_transition.to_csv(f"{OUTDIR}/cluster_transition_matrix_{analysis_type}.csv")
print(f"Saved: {OUTDIR}/cluster_transition_matrix_{analysis_type}.csv")

# ============================================================
# 9. まとめレポート
# ============================================================
print(f"\n{'='*60}")
print("CYCLE DETECTION SUMMARY")
print(f"{'='*60}")
print(f"Analysis type: {analysis_type}")
print(f"Number of clusters: {n_clusters}")
print(f"Graph edges: {G.number_of_edges()}")
print(f"Detected cycles: {len(cycles)}")
print(f"Strongly connected components: {len(scc)}")

if len(cycle_df) > 0:
    print("\nCycle details:")
    for _, row in cycle_df.iterrows():
        print(f"  Cycle {int(row['cycle_id'])}: {row['states']} "
              f"(length: {int(row['length'])}, prob: {row['cycle_probability']:.4f})")

print(f"\nOutput directory: {OUTDIR}/")
print("Files generated:")
print(f"  - cluster_transition_cycles_{analysis_type}.png")
print(f"  - detected_cycles_{analysis_type}.csv")
print(f"  - strongly_connected_components_{analysis_type}.csv")
print(f"  - cluster_transition_matrix_{analysis_type}.csv")

print(f"\n{'='*60}")

# ============================================================
# 10. 生物学的解釈
# ============================================================
if len(cycles) > 0:
    print("\nBIOLOGICAL INTERPRETATION:")
    print("- Cycles indicate potential differentiation loops or feedback mechanisms")
    print("- High cycle probability suggests strong bidirectional transitions")
    print("- Multiple cycles may indicate complex cell fate decisions")
    print("- Consider batch effects if unexpected cycles are detected")

print("\nAnalysis complete!")