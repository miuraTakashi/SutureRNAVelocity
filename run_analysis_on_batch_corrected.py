"""
バッチ補正済みデータを使ってvelocity分析を行うスクリプト

- integrated_velocity_batch_corrected.h5ad を読み込み
- OGクラスタをフィルタ
- Velocityグラフ計算
- ストリームプロット
- サイクル検出（velocityのみ）
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
import networkx as nx

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

adata_path = "integrated_velocity_batch_corrected_fixed.h5ad"
if not os.path.exists(adata_path):
    raise FileNotFoundError(f"{adata_path} not found.")

print(f"Loading {adata_path}...")
adata = sc.read_h5ad(adata_path)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# OG/PO クラスタラベルの整備（仮定: og_cluster 列がある）
if "og_cluster" not in adata.obs:
    print("Warning: og_cluster not found. Assuming all are OG clusters.")
    adata.obs["og_cluster"] = "OG"  # 仮定
else:
    adata.obs["og_cluster"] = pd.Categorical(
        adata.obs["og_cluster"].astype(str),
        categories=["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"],
    )

clusters = ["OG1", "OG2", "OG3"]
if "og_cluster" in adata.obs and all(c in adata.obs["og_cluster"].unique() for c in clusters):
    adata = adata[adata.obs["og_cluster"].isin(clusters)].copy()
    print(f"Filtering to clusters: {clusters}")
    print(f"  {adata.n_obs} cells remain after filtering.")
else:
    print("OG clusters not found. Proceeding with all data.")
    clusters = adata.obs["og_cluster"].unique() if "og_cluster" in adata.obs else ["all"]

# データの型を修正
print("Converting layers to float...")
for key in adata.layers:
    adata.layers[key] = adata.layers[key].astype(np.float32)

# Velocityグラフ計算
if "velocity_graph" not in adata.uns:
    print("Computing velocity graph...")
    try:
        # Use dynamical model like in run_og_exact.py
        scv.tl.recover_dynamics(adata, n_jobs=4)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
    except Exception as e:
        print(f"Error computing velocity graph: {e}")
        print("Skipping velocity analysis.")
        exit(1)
else:
    print("Velocity graph already exists.")

# UMAP再計算（バッチ補正後）
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

# ストリームプロット
outdir = "figures_batch_corrected_analysis"
os.makedirs(outdir, exist_ok=True)
sc.settings.figdir = outdir

OG_COLORS = {
    "OG1": "#e41a1c",
    "OG2": "#ff7f00",
    "OG3": "#4daf4a",
    "OG4": "#984ea3",
    "PO1": "#377eb8",
    "PO2": "#a65628",
}
OG_PALETTE = list(OG_COLORS.values())

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="og_cluster",
    palette=OG_PALETTE,
    title="Velocity Stream (Batch Corrected Integrated)",
    save="_velocity_stream_batch_corrected_integrated.png",
)

# サイクル検出（velocityのみ）
print("Detecting cycles using velocity only...")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
T = vk.transition_matrix

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

print("Cluster transition matrix (velocity only):")
print(cluster_flux.round(4).to_string())

threshold = cluster_flux.values.mean() * 0.5
print(f"Using threshold {threshold:.4f} for cycle detection.")
G = nx.DiGraph()
for a in clusters:
    for b in clusters:
        if a == b:
            continue
        weight = cluster_flux.loc[a, b]
        if weight > threshold:
            G.add_edge(a, b, weight=float(weight))

cycles = list(nx.simple_cycles(G))
if len(cycles) == 0:
    print("No cycles detected.")
else:
    print(f"Detected {len(cycles)} cycle(s):")
    for cycle in cycles:
        print("  -> ".join(cycle) + " -> " + cycle[0])

# 可視化
plt.figure(figsize=(8, 8))
pos = nx.circular_layout(G)
edge_weights = [data["weight"] for _, _, data in G.edges(data=True)]
max_weight = max(edge_weights) if edge_weights else 1.0

nx.draw_networkx_nodes(G, pos, node_color="#1f77b4", node_size=1100)
nx.draw_networkx_labels(G, pos, font_size=10, font_color="white")

all_edges = list(G.edges(data=True))
edge_colors = ["lightgray"] * len(all_edges)
edge_widths = [1.0 + 4.0 * (data["weight"] / max_weight) for _, _, data in all_edges]

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

edge_labels = {(u, v): f"{data['weight']:.3f}" for u, v, data in all_edges}
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

plt.title("Cluster transition cycles (Velocity only, Batch Corrected)")
plt.axis("off")
plot_path = os.path.join(outdir, "cluster_cycle_strengths_batch_corrected.png")
plt.tight_layout()
plt.savefig(plot_path, dpi=200)
plt.close()
print(f"Saved cycle plot to {plot_path}")

print("Done.")