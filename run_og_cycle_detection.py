"""
OG/PO クラスタ間の遷移サイクルを検出するスクリプト

- E17_og_exact_velocity.h5ad を読み込み
- CellRank の VelocityKernel + ConnectivityKernel で遷移確率行列を構成
- OG/PO クラスタ間の正味遷移を集計
- NetworkX でクラスタレベルの有向サイクルを検出
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import cellrank as cr
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

adata_path = "E17_og_exact_velocity.h5ad"
if not os.path.exists(adata_path):
    raise FileNotFoundError(f"{adata_path} not found. Please run the preprocessing pipeline first.")

print(f"Loading {adata_path}...")
adata = sc.read_h5ad(adata_path)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# OG/PO クラスタラベルの整備 （OG1-3 のみを対象）
adata.obs["og_cluster"] = pd.Categorical(
    adata.obs["og_cluster"].astype(str),
    categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
)
clusters = ["OG1", "OG2", "OG3"]
adata = adata[adata.obs["og_cluster"].isin(clusters)].copy()
print(f"Filtering to clusters: {clusters}")
print(f"  {adata.n_obs} cells remain after filtering.")

# velocity_graph が存在しない場合は再計算
if "velocity_graph" not in adata.uns:
    print("Computing velocity graph...")
    scv.tl.velocity_graph(adata)

print("Building CellRank transition kernels...")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
T = combined_kernel.transition_matrix

# クラスタ対クラスタの平均遷移確率を計算
print("Computing cluster-level transition flux...")
cluster_flux = pd.DataFrame(0.0, index=clusters, columns=clusters)
for a in clusters:
    idx_a = np.where(adata.obs["og_cluster"].values == a)[0]
    if idx_a.size == 0:
        continue
    for b in clusters:
        idx_b = np.where(adata.obs["og_cluster"].values == b)[0]
        if idx_b.size == 0:
            continue
        sub = T[idx_a, :][:, idx_b]
        cluster_flux.loc[a, b] = np.asarray(sub.sum(axis=1)).mean()

print("Cluster transition matrix (mean outgoing probability):")
print(cluster_flux.round(4).to_string())

# 有向グラフを作成し、サイクルを検出
threshold = cluster_flux.values.mean() * 0.5
print(f"Using threshold {threshold:.4f} to include edges in the cycle graph.")
G = nx.DiGraph()
for a in clusters:
    for b in clusters:
        if a == b:
            continue
        weight = cluster_flux.loc[a, b]
        if weight > threshold:
            G.add_edge(a, b, weight=float(weight))

print("Edges used for cycle detection:")
for u, v, data in G.edges(data=True):
    print(f"  {u} -> {v} (weight={data['weight']:.4f})")

cycles = list(nx.simple_cycles(G))
if len(cycles) == 0:
    print("No directed cluster cycles detected in the thresholded graph.")
else:
    print(f"Detected {len(cycles)} directed cluster cycle(s):")
    for cycle in cycles:
        print("  -> ".join(cycle) + " -> " + cycle[0])

# 強連結成分も確認
scc = [c for c in nx.strongly_connected_components(G) if len(c) > 1]
if len(scc) > 0:
    print("\nStrongly connected cluster components (potential cycle modules):")
    for comp in scc:
        print("  " + ", ".join(sorted(comp)))
else:
    print("\nNo multi-node strongly connected components found.")

# サイクルの可視化
outdir = "cycle_detection"
os.makedirs(outdir, exist_ok=True)

plt.figure(figsize=(8, 8))
pos = nx.circular_layout(G)
edge_weights = [data["weight"] for _, _, data in G.edges(data=True)]
max_weight = max(edge_weights) if edge_weights else 1.0

# 全エッジを描画
nx.draw_networkx_nodes(G, pos, node_color="#1f77b4", node_size=1100)
nx.draw_networkx_labels(G, pos, font_size=10, font_color="white")

all_edges = list(G.edges(data=True))
edge_colors = ["lightgray"] * len(all_edges)
edge_widths = [1.0 + 4.0 * (data["weight"] / max_weight) for _, _, data in all_edges]

# サイクル内エッジは赤で強調
cycle_edges = set()
for cycle in cycles:
    for i in range(len(cycle)):
        a = cycle[i]
        b = cycle[(i + 1) % len(cycle)]
        cycle_edges.add((a, b))
for idx, (u, v, data) in enumerate(all_edges):
    if (u, v) in cycle_edges:
        edge_colors[idx] = "red"
        edge_widths[idx] = edge_widths[idx] + 2.0

nx.draw_networkx_edges(
    G,
    pos,
    edgelist=[(u, v) for u, v, _ in all_edges],
    width=edge_widths,
    edge_color=edge_colors,
    arrows=True,
    arrowstyle='-|>',
    arrowsize=20,
    connectionstyle="arc3,rad=0.1",
)

# エッジ重みの注釈
edge_labels = {(u, v): f"{data['weight']:.3f}" for u, v, data in all_edges}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

plt.title("OG/PO cluster transition cycles and edge strengths")
plt.axis("off")
plot_path = os.path.join(outdir, "cluster_cycle_strengths.png")
plt.tight_layout()
plt.savefig(plot_path, dpi=200)
plt.close()
print(f"Saved cycle visualization to {plot_path}")

print("Done.")
