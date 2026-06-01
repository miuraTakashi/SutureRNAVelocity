"""
OG/POクラスタ図の上に疑似時系列を重ねて表示するスクリプト
"""
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR = Path(__file__).resolve().parent
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

cluster_colors = {
    "OG1": "#FF0000",
    "OG2": "#FF8000",
    "OG3": "#FFFF00",
    "OG4": "#00FF00",
    "PO1": "#0080FF",
    "PO2": "#8000FF"
}

h5ad_files = [
    "E17_batch_1_possorted_genome_bam_D5AOD_velocity_OG_filtered_v2.h5ad",
    "E17_batch_2_possorted_genome_bam_Z9D92_velocity_OG_filtered_v2.h5ad",
    "E17_batch_3_possorted_genome_bam_WACA8_velocity_OG_filtered_v2.h5ad"
]

for h5ad_path in h5ad_files:
    if not Path(h5ad_path).exists():
        print(f"File not found: {h5ad_path}")
        continue

    sample_name = Path(h5ad_path).stem.replace("_velocity_OG_filtered_v2", "")
    print(f"\n{'='*60}")
    print(f"Processing {sample_name}")
    print('='*60)

    adata = sc.read_h5ad(h5ad_path)
    if "cluster" not in adata.obs:
        print("  No cluster info found. Skipping.")
        continue

    if "X_umap" not in adata.obsm:
        print("  Recomputing UMAP...")
        if "X_pca" not in adata.obsm:
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))
        sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
        sc.tl.umap(adata)
        print("  UMAP recomputed.")

    if "dpt_pseudotime" not in adata.obs:
        print("  Computing pseudotime (DPT)...")
        try:
            if "X_pca" not in adata.obsm:
                sc.pp.scale(adata, max_value=10)
                sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))
            sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
            sc.tl.diffmap(adata)

            if "cluster" in adata.obs and "OG1" in adata.obs["cluster"].unique():
                root_name = adata.obs[adata.obs["cluster"] == "OG1"].index[0]
            else:
                root_name = adata.obs_names[0]
            adata.uns["iroot"] = adata.obs_names.get_loc(root_name)

            sc.tl.dpt(adata)
            if "dpt_pseudotime" in adata.obs:
                print("  DPT pseudotime computed.")
            else:
                raise RuntimeError("DPT did not populate dpt_pseudotime")
        except Exception as e:
            print(f"  DPT failed: {e}")
            if "velocity_pseudotime" in adata.obs:
                adata.obs["dpt_pseudotime"] = adata.obs["velocity_pseudotime"]
                print("  velocity_pseudotime copied to dpt_pseudotime.")
            elif "latent_time" in adata.obs:
                adata.obs["dpt_pseudotime"] = adata.obs["latent_time"]
                print("  latent_time copied to dpt_pseudotime.")
            else:
                try:
                    if "neighbors" not in adata.uns or adata.obsp.get('distances') is None:
                        if "X_pca" not in adata.obsm:
                            sc.pp.scale(adata, max_value=10)
                            sc.tl.pca(adata, n_comps=min(30, adata.shape[1]-1))
                        sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm['X_pca'].shape[1]))
                    scv.tl.velocity_pseudotime(adata)
                    adata.obs["dpt_pseudotime"] = adata.obs["velocity_pseudotime"]
                    print("  Velocity pseudotime computed.")
                except Exception as e2:
                    print(f"  velocity_pseudotime failed: {e2}")
                    continue

    pseudotime = adata.obs["dpt_pseudotime"].astype(float)
    umap = adata.obsm["X_umap"]

    fig, ax = plt.subplots(figsize=(12, 10))
    unique_clusters = sorted(adata.obs["cluster"].unique(), key=lambda x: str(x))

    # Base cluster plot
    for cluster in unique_clusters:
        mask = adata.obs["cluster"] == cluster
        color = cluster_colors.get(cluster, "#999999")
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=25,
            c=color,
            alpha=0.35,
            label=cluster,
            edgecolors="none"
        )

    # Pseudotime overlay
    scatter = ax.scatter(
        umap[:, 0],
        umap[:, 1],
        c=pseudotime,
        cmap="viridis",
        s=35,
        alpha=0.9,
        edgecolors="none"
    )

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.7)
    cbar.set_label("pseudotime", fontsize=12)

    ax.set_title(f"{sample_name} - Cluster + Pseudotime", fontsize=16, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(False)

    out_file = OUT_DIR / f"{sample_name}_cluster_pseudotime_overlay.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved overlay plot: {out_file.name}")

print("\nAll pseudotime overlay plots generated.")
