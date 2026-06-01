"""
サイクル検出デモンストレーション
人工的な遷移確率行列を使ってサイクル検出アルゴリズムを示す
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

# ============================================================
# 1. 人工的な遷移確率行列の作成（OG1-4 + PO1-2）
# ============================================================
clusters = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]

print("Creating artificial transition matrix with cycles...")

# 基本的な分化経路 + サイクルを含む行列
transition_matrix = pd.DataFrame(0.0, index=clusters, columns=clusters)

# 基本分化経路: OG1 → OG2 → OG3 → OG4
transition_matrix.loc["OG1", "OG2"] = 0.6
transition_matrix.loc["OG2", "OG3"] = 0.5
transition_matrix.loc["OG3", "OG4"] = 0.4

# サイクル1: OG4 → OG1 (脱分化)
transition_matrix.loc["OG4", "OG1"] = 0.3

# サイクル2: OG3 ↔ PO1 (分岐)
transition_matrix.loc["OG3", "PO1"] = 0.2
transition_matrix.loc["PO1", "OG3"] = 0.4

# PO系列: PO1 → PO2
transition_matrix.loc["PO1", "PO2"] = 0.3

# 自己遷移（各状態に留まる確率）
for cl in clusters:
    # 他の遷移確率の合計を計算
    total_out = transition_matrix.loc[cl].sum()
    # 残りを自己遷移に割り当て
    transition_matrix.loc[cl, cl] = max(0.1, 1.0 - total_out)

# 正規化（各行の合計が1になるように）
for cl in clusters:
    row_sum = transition_matrix.loc[cl].sum()
    if row_sum > 0:
        transition_matrix.loc[cl] = transition_matrix.loc[cl] / row_sum

print("Artificial transition matrix:")
print(transition_matrix.round(3))

# ============================================================
# 2. NetworkX グラフの作成
# ============================================================
print("\nCreating NetworkX graph...")

G = nx.DiGraph()

# ノード追加
og_colors = {
    "OG1": "#e41a1c", "OG2": "#ff7f00", "OG3": "#4daf4a",
    "OG4": "#984ea3", "PO1": "#377eb8", "PO2": "#a65628"
}

for cl in clusters:
    G.add_node(cl, color=og_colors[cl])

# エッジ追加（遷移確率 > 0.05 の場合）
threshold = 0.05
for source in clusters:
    for target in clusters:
        prob = transition_matrix.loc[source, target]
        if prob > threshold:
            G.add_edge(source, target, weight=prob, probability=prob)

print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ============================================================
# 3. サイクル検出
# ============================================================
print("\nDetecting cycles...")

cycles = list(nx.simple_cycles(G))
print(f"Found {len(cycles)} cycles")

cycle_info = []
for i, cycle in enumerate(cycles):
    cycle_length = len(cycle)

    # サイクル内の遷移確率の積
    cycle_prob = 1.0
    for j in range(cycle_length):
        source = cycle[j]
        target = cycle[(j + 1) % cycle_length]
        prob = transition_matrix.loc[source, target]
        cycle_prob *= prob

    cycle_info.append({
        'cycle_id': i + 1,
        'cycle': cycle,
        'length': cycle_length,
        'cycle_probability': cycle_prob,
        'states': ' → '.join(cycle)
    })

cycle_df = pd.DataFrame(cycle_info)
print("\nDetected cycles:")
if len(cycle_df) > 0:
    print(cycle_df.to_string(index=False))
else:
    print("No cycles detected.")

# ============================================================
# 4. 強連結成分分析
# ============================================================
print("\nAnalyzing strongly connected components...")

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
                prob = transition_matrix.loc[source, target]
                internal_probs.append(prob)

    avg_internal_prob = np.mean(internal_probs) if internal_probs else 0.0

    scc_info.append({
        'scc_id': i + 1,
        'size': component_size,
        'states': ', '.join(component_list),
        'avg_internal_probability': avg_internal_prob,
        'is_cycle': component_size > 1
    })

scc_df = pd.DataFrame(scc_info)
print("\nStrongly connected components:")
print(scc_df.to_string(index=False))

# ============================================================
# 5. 可視化
# ============================================================
print("\nCreating visualization...")

plt.figure(figsize=(12, 8))

# ノード位置（spring layout）
pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

# ノード描画
node_colors = [og_colors[node] for node in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=2000, alpha=0.8)

# エッジ描画
edges = G.edges()
edge_weights = [G[u][v]['weight'] * 10 for u, v in edges]
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

    nx.draw_networkx_edges(G, pos, edgelist=cycle_edges, width=4,
                          edge_color='red', alpha=0.9, arrows=True,
                          arrowsize=25, connectionstyle='arc3,rad=0.2')

# エッジ確率ラベル
edge_labels = {}
for u, v in G.edges():
    prob = G[u][v]['probability']
    if prob > threshold:
        edge_labels[(u, v)] = f'{prob:.2f}'

nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, alpha=0.8)

plt.title('Artificial Transition Graph with Detected Cycles', fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.savefig('artificial_cycles_demo.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: artificial_cycles_demo.png")

# ============================================================
# 6. CellRank GPCCA シミュレーション
# ============================================================
print("\nSimulating CellRank GPCCA analysis...")

# 遷移確率行列をnumpy配列に変換
P = transition_matrix.values

# 固有値分解（簡易版）
eigenvals, eigenvecs = np.linalg.eig(P.T)  # 左固有ベクトル用に転置

# 最大固有値（1に最も近い）に対応する固有ベクトル
max_eigenval_idx = np.argmin(np.abs(eigenvals - 1.0))
stationary_dist = np.real(eigenvecs[:, max_eigenval_idx])
stationary_dist = stationary_dist / stationary_dist.sum()  # 正規化

print("Stationary distribution (π):")
for i, cl in enumerate(clusters):
    print(".3f")

# Terminal state の候補（定常分布が大きい順）
terminal_candidates = sorted(zip(clusters, stationary_dist),
                           key=lambda x: x[1], reverse=True)

print("\nTerminal state candidates (by stationary distribution):")
for i, (state, prob) in enumerate(terminal_candidates[:3]):
    print(".3f")

# ============================================================
# 7. まとめ
# ============================================================
print(f"\n{'='*60}")
print("CYCLE DETECTION DEMONSTRATION SUMMARY")
print(f"{'='*60}")
print(f"Clusters: {len(clusters)}")
print(f"Graph edges: {G.number_of_edges()}")
print(f"Detected cycles: {len(cycles)}")
print(f"Strongly connected components: {len(scc)}")

if len(cycles) > 0:
    print("\nCycle analysis:")
    for _, row in cycle_df.iterrows():
        print(f"  Cycle {int(row['cycle_id'])}: {row['states']}")
        print(".4f")

print("\nKey findings:")
print("- Cycles represent potential feedback loops in cell differentiation")
print("- Stationary distribution helps identify terminal states within cycles")
print("- GPCCA would cluster cyclic components and select most stable states")
print("- In real data, cycles might indicate dedifferentiation or batch effects")

print(f"\nOutput: artificial_cycles_demo.png")
print(f"{'='*60}")