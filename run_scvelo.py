"""
RNA Velocity analysis of mouse cranial suture E17 scRNA-seq (GSE163693)
Using scVelo with velocyto-generated loom files.
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import scvelo as scv
import anndata as ad
import loompy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_loom_as_anndata(filename):
    """Read a velocyto loom file into AnnData (loompy backend)."""
    with loompy.connect(filename, "r") as ds:
        X = sp.csc_matrix(ds[:, :]).T  # cells x genes
        layers = {}
        for key in ds.layers.keys():
            if key == "":
                continue
            layers[key] = sp.csc_matrix(ds.layers[key][:, :]).T
        obs_names = pd.Index(ds.col_attrs["CellID"][:])
        var_names = pd.Index(ds.row_attrs["Gene"][:])
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=var_names),
        )
        for key, mat in layers.items():
            adata.layers[key] = mat
    return adata

warnings.filterwarnings("ignore")
scv.settings.verbosity = 3
scv.settings.presenter_view = True
sc.settings.figdir = "figures"
os.makedirs("figures", exist_ok=True)

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

# ============================================================
# 1. Load GEO expression matrix (filtered cells)
# ============================================================
print("=" * 60)
print("Step 1: Loading GEO expression matrix")
print("=" * 60)

adata = sc.read_10x_mtx(
    "GSE163693_RAW/",
    prefix="GSM4983999_E17_",
    var_names="gene_symbols",
    cache=True,
)
adata.var_names_make_unique()
print(f"  GEO matrix: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ============================================================
# 2. Standard preprocessing (Scanpy)
# ============================================================
print("=" * 60)
print("Step 2: Preprocessing")
print("=" * 60)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
print(f"  After QC: {adata.shape[0]} cells x {adata.shape[1]} genes")

adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)

sc.pl.umap(adata, color=["leiden"], save="_clusters.png")
print(f"  Found {adata.obs['leiden'].nunique()} clusters")

# ============================================================
# 3. Load velocyto loom files and merge
# ============================================================
print("=" * 60)
print("Step 3: Loading loom files")
print("=" * 60)

loom_files = sorted(glob.glob("loom_output/*.loom"))
if not loom_files:
    raise FileNotFoundError(
        "No loom files found in loom_output/. Run velocyto first (run_velocyto.sh)."
    )

loom_list = []
for f in loom_files:
    print(f"  Loading {f}")
    ldata = read_loom_as_anndata(f)
    ldata.var_names_make_unique()
    loom_list.append(ldata)

if len(loom_list) == 1:
    ldata_merged = loom_list[0]
else:
    ldata_merged = ad.concat(loom_list, join="outer")

print(f"  Loom data: {ldata_merged.shape[0]} cells x {ldata_merged.shape[1]} genes")

# ============================================================
# 4. Merge velocity data into expression AnnData
# ============================================================
print("=" * 60)
print("Step 4: Merging velocity data")
print("=" * 60)

# Clean barcode names for matching
# velocyto barcodes: "samplename:BARCODE" or just "BARCODE"
# GEO barcodes: "AAACCTGAGACAATAC-1"

loom_barcodes = ldata_merged.obs_names.tolist()
if ":" in loom_barcodes[0]:
    ldata_merged.obs_names = [bc.split(":")[-1] for bc in loom_barcodes]
    # Remove any trailing "x" added by velocyto
    ldata_merged.obs_names = [bc.replace("x", "") if bc.endswith("x") else bc for bc in ldata_merged.obs_names]

ldata_merged.obs_names_make_unique()
ldata_merged.var_names_make_unique()

common_cells = adata.obs_names.intersection(ldata_merged.obs_names)
print(f"  Common cells: {len(common_cells)} / {adata.shape[0]} (GEO) vs {ldata_merged.shape[0]} (loom)")

if len(common_cells) == 0:
    print("  WARNING: No matching barcodes. Trying to strip suffix...")
    adata_bc_base = [bc.rsplit("-", 1)[0] for bc in adata.obs_names]
    loom_bc_base = [bc.rsplit("-", 1)[0] for bc in ldata_merged.obs_names]

    adata.obs["barcode_base"] = adata_bc_base
    ldata_merged.obs["barcode_base"] = loom_bc_base

    common_base = set(adata_bc_base) & set(loom_bc_base)
    print(f"  Common base barcodes: {len(common_base)}")

adata_vel = scv.utils.merge(adata, ldata_merged)
print(f"  Merged: {adata_vel.shape[0]} cells x {adata_vel.shape[1]} genes")

# ============================================================
# 5. scVelo RNA velocity (stochastic model)
# ============================================================
print("=" * 60)
print("Step 5: RNA velocity estimation (stochastic)")
print("=" * 60)

scv.pp.filter_and_normalize(adata_vel, min_shared_counts=20)
sc.pp.highly_variable_genes(adata_vel, n_top_genes=2000, flavor="seurat")
# Recompute neighbors in scvelo format to avoid graph format mismatch
scv.pp.moments(adata_vel, n_pcs=30, n_neighbors=30)

print("  Skipping stochastic model (known numpy compatibility issue).")
print("  Proceeding directly to dynamical model.")

# ============================================================
# 6. Dynamical model (higher accuracy)
# ============================================================
print("=" * 60)
print("Step 6: RNA velocity estimation (dynamical)")
print("=" * 60)

scv.tl.recover_dynamics(adata_vel, n_jobs=4)
scv.tl.velocity(adata_vel, mode="dynamical")
scv.tl.velocity_graph(adata_vel)

scv.pl.velocity_embedding_stream(
    adata_vel, basis="umap", color="leiden",
    save="velocity_stream_dynamical.png",
)

# ============================================================
# 7. Latent time and driver genes
# ============================================================
print("=" * 60)
print("Step 7: Latent time and driver genes")
print("=" * 60)

scv.tl.latent_time(adata_vel)
scv.pl.scatter(
    adata_vel, color="latent_time", color_map="gnuplot",
    save="latent_time.png",
)

top_genes = adata_vel.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_vel, var_names=top_genes[:50], sortby="latent_time", save="heatmap_top_genes.png")

scv.tl.rank_velocity_genes(adata_vel, groupby="leiden", min_corr=0.3)
df_velocity_genes = pd.DataFrame(adata_vel.uns["rank_velocity_genes"]["names"])
df_velocity_genes.to_csv("velocity_driver_genes.csv", index=False)
print("  Driver genes saved to velocity_driver_genes.csv")
print(df_velocity_genes.head(10))

# ============================================================
# 8. Save results
# ============================================================
print("=" * 60)
print("Step 8: Saving results")
print("=" * 60)

adata_vel.write("E17_suture_velocity.h5ad")
print("  Saved to E17_suture_velocity.h5ad")

scv.pl.proportions(adata_vel, groupby="leiden", save="spliced_proportions.png")

print("=" * 60)
print("DONE! All results saved.")
print(f"  - Figures: figures/")
print(f"  - AnnData: E17_suture_velocity.h5ad")
print(f"  - Driver genes: velocity_driver_genes.csv")
print("=" * 60)
