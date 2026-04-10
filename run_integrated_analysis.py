"""
総合分析: ストリームプロット + サイクル検出 + CellRank
バッチ補正済みvelocityデータを使用
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

# ============================================================
# 1. Load data
# ============================================================
print("Loading data...")
adata = sc.read_h5ad("E17_og_integrated_harmony_filtered_velocity.h5ad")
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Clusters: {adata.obs['og_cluster'].value_counts().to_dict()}")

outdir = "figures_integrated_analysis"
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

# ============================================================
# 2. Stream plot
# ============================================================
print("\n=== Generating Stream Plot ===")
try:
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="og_cluster",
        palette=OG_COLORS,
        title="RNA Velocity Stream (Batch Corrected, OG1-4/PO1-2)",
        save="_velocity_stream.png",
        dpi=150,
    )
    print("Stream plot saved.")
except Exception as e:
    print(f"Error: {e}")

# ============================================================
# 3. Cycle detection (velocity-based)
# ============================================================
print("\n=== Cycle Detection ===")
try:
    clusters = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
    clusters = [c for c in clusters if c in adata.obs["og_cluster"].unique()]
    
    # Build transition matrix from velocity
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    T = vk.transition_matrix
    
    # Cluster-level transition flux
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
    
    print("Cluster transition matrix:")
    print(cluster_flux.round(4).to_string())
    
    # Cycle detection
    threshold = cluster_flux.values.mean() * 0.5
    G = nx.DiGraph()
    for a in clusters:
        for b in clusters:
            if a == b:
                continue
            weight = cluster_flux.loc[a, b]
            if weight > threshold:
                G.add_edge(a, b, weight=float(weight))
    
    cycles = list(nx.simple_cycles(G))
    print(f"\nDetected {len(cycles)} cycle(s):")
    for cycle in cycles:
        print("  " + " -> ".join(cycle) + " -> " + cycle[0])
    
    # Visualize cycles
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
        G, pos,
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
    
    plt.title("Cluster transition cycles (velocity-based)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"{outdir}/cycle_detection.png", dpi=200)
    plt.close()
    print(f"Cycle plot saved to {outdir}/cycle_detection.png")
    
except Exception as e:
    print(f"Error in cycle detection: {e}")

# ============================================================
# 4. CellRank analysis
# ============================================================
print("\n=== CellRank Analysis ===")
try:
    # Build combined kernel
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    combined_kernel = 0.8 * vk + 0.2 * ck
    
    # Compute macrostates with GPCCA
    print("Computing macrostates...")
    combined_kernel.compute_macrostates(n_states=4)
    print(f"Macrostates: {adata.obs['macrostates'].unique()}")
    
    # Compute terminal states
    print("Computing terminal states...")
    combined_kernel.compute_terminal_states()
    print(f"Terminal states: {adata.obs['terminal_states'].unique()}")
    
    # Terminal state probabilities
    print("Computing terminal state probabilities...")
    combined_kernel.compute_absorption_probabilities()
    
    print("CellRank analysis completed.")
    
    # Visualization
    scv.pl.scatter(
        adata,
        basis="umap",
        color="macrostates",
        title="CellRank Macrostates",
        save="_macrostates.png",
    )
    
    scv.pl.scatter(
        adata,
        basis="umap",
        color="terminal_states",
        title="CellRank Terminal States",
        save="_terminal_states.png",
    )
    
    print("CellRank plots saved.")
    
except Exception as e:
    print(f"Error in CellRank: {e}")

# ============================================================
# 5. Save results
# ============================================================
adata.write(f"{outdir}/analysis_results.h5ad")
print(f"\nAnalysis completed. Results saved to {outdir}/")
