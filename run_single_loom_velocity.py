"""
3つのloomファイルを個別に RNA velocity 解析し、結果を保存するスクリプト。
"""
import os
from pathlib import Path
import scanpy as sc
import scvelo as scv
import numpy as np
import warnings

warnings.filterwarnings("ignore")

WORKDIR = Path(__file__).resolve().parent
LOOM_DIR = Path("/home/user/share/SutureRNAVelocity/loom_output")
OUT_DIR = WORKDIR / "figures_single_loom_velocity"
OUT_DIR.mkdir(exist_ok=True)

loom_files = sorted(LOOM_DIR.glob("*.loom"))
if not loom_files:
    raise FileNotFoundError(f"No loom files found in {LOOM_DIR}")

for loom_path in loom_files:
    sample_name = loom_path.stem
    print(f"\nProcessing {sample_name}")

    adata = sc.read_loom(str(loom_path), sparse=False)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Ensure layers exist and are numeric
    for layer in ["spliced", "unspliced", "ambiguous"]:
        if layer in adata.layers:
            adata.layers[layer] = adata.layers[layer].astype(np.float32)
    if adata.X is not None:
        adata.X = adata.X.astype(np.float32)

    # Basic preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    hvgs = adata.var.highly_variable
    if hvgs.sum() < 100:
        raise ValueError(f"Too few HVGs for {sample_name}: {hvgs.sum()}")
    adata = adata[:, hvgs].copy()

    # Save a copy before velocity processing to preserve UMAP
    print("Computing UMAP...")
    try:
        adata_for_umap = adata.copy()
        sc.pp.scale(adata_for_umap, max_value=10)
        sc.tl.pca(adata_for_umap, n_comps=min(30, adata_for_umap.shape[1]-1))
        sc.pp.neighbors(adata_for_umap, n_pcs=min(30, adata_for_umap.obsm['X_pca'].shape[1]))
        sc.tl.umap(adata_for_umap)
        
        # Save UMAP to main object
        adata.obsm['umap'] = adata_for_umap.obsm['umap']
        if 'X_pca' in adata_for_umap.obsm:
            adata.obsm['X_pca'] = adata_for_umap.obsm['X_pca']
        if 'distances' in adata_for_umap.obsp:
            adata.obsp['distances'] = adata_for_umap.obsp['distances']
        if 'connectivities' in adata_for_umap.obsp:
            adata.obsp['connectivities'] = adata_for_umap.obsp['connectivities']
        if 'umap' in adata_for_umap.uns:
            adata.uns['umap'] = adata_for_umap.uns['umap']
        if 'neighbors' in adata_for_umap.uns:
            adata.uns['neighbors'] = adata_for_umap.uns['neighbors']
        print("  UMAP computed and saved")
    except Exception as e:
        print(f"  UMAP computation failed: {e}")

    # Velocity processing
    print("Computing velocity...")
    try:
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        scv.tl.velocity(adata, mode="stochastic")
        scv.tl.velocity_graph(adata)
        print("  Velocity computed successfully")
    except Exception as e:
        print(f"  Stochastic velocity failed: {e}")
        try:
            print("  Trying dynamical model...")
            adata_dyn = adata.copy()
            scv.tl.recover_dynamics(adata_dyn, n_jobs=4)
            scv.tl.velocity(adata_dyn, mode="dynamical")
            scv.tl.velocity_graph(adata_dyn)
            adata = adata_dyn
            print("  Dynamical model succeeded")
        except Exception as e2:
            print(f"  Dynamical model failed: {e2}")
            print("  Skipping velocity graph for this sample")

    # Save AnnData and plots
    out_h5ad = WORKDIR / f"{sample_name}_velocity.h5ad"
    adata.write(out_h5ad)
    print(f"  Saved AnnData: {out_h5ad}")

    try:
        ax = scv.pl.velocity_embedding_stream(adata, basis="umap", color=None, show=False, legend_loc="none")
        ax.figure.savefig(OUT_DIR / f"{sample_name}_velocity_stream.png", dpi=300, bbox_inches="tight")
        ax.figure.clf()
        print(f"  Saved velocity stream plot: {OUT_DIR / f'{sample_name}_velocity_stream.png'}")
    except Exception as e:
        print(f"  Velocity stream plot failed: {e}")

    import matplotlib.pyplot as plt
    sc.pl.umap(adata, color=None, show=False)
    plt.gcf().savefig(OUT_DIR / f"{sample_name}_umap.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved UMAP plot: {OUT_DIR / f'{sample_name}_umap.png'}")

print("\nAll loom files processed.")
