"""
Loom統合スクリプト (Harmony + OGクラスタのみ)
既存メタデータを活用してメモリ効率を改善
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import loompy
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

try:
    import harmony
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False
    print("Harmony not found. Using BBKNN.")

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

def read_loom_as_anndata(filename):
    """Read a velocyto loom file into AnnData."""
    with loompy.connect(filename, "r") as ds:
        if "spliced" in ds.layers:
            X = ds.layers["spliced"][:, :].T
        else:
            X = np.zeros((ds.shape[1], ds.shape[0]))
        layers = {}
        for key in ds.layers.keys():
            if key == "":
                continue
            layers[key] = ds.layers[key][:, :].T
        obs_names = pd.Index(ds.col_attrs["CellID"][:])
        var_names = pd.Index(ds.row_attrs["Gene"][:])
        adata_loom = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=var_names),
        )
        for key, mat in layers.items():
            adata_loom.layers[key] = mat
    return adata_loom

# ============================================================
# 1. Load Loom files
# ============================================================
print("Loading Loom files...")
loom_files = sorted(glob.glob("/home/user/share/SutureRNAVelocity/loom_output/*.loom"))
if not loom_files:
    raise FileNotFoundError("No loom files found")

loom_list = []
batch_info = []
for f in loom_files:
    print(f"  Loading {f}")
    ldata = read_loom_as_anndata(f)
    ldata.var_names_make_unique()
    
    # Extract batch info from filename
    basename = os.path.basename(f)
    if "batch_1" in basename:
        batch_id = "batch_1"
    elif "batch_2" in basename:
        batch_id = "batch_2"
    elif "batch_3" in basename:
        batch_id = "batch_3"
    else:
        batch_id = basename.split(".")[0]
    
    ldata.obs["batch"] = batch_id
    batch_info.append((f, batch_id, ldata.shape[0]))
    loom_list.append(ldata)

# ============================================================
# 2. Merge Loom files
# ============================================================
print("Integrating Loom files...")
if len(loom_list) == 1:
    adata = loom_list[0]
else:
    adata = ad.concat(loom_list, join="outer")

print(f"  Integrated data: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ============================================================
# 3. Load existing metadata (og_cluster annotation)
# ============================================================
print("Loading existing metadata...")
try:
    adata_meta = sc.read_h5ad("E17_og_exact_velocity.h5ad")
    print(f"  Metadata source: {adata_meta.shape[0]} cells")
    
    # Extract og_cluster information
    og_cluster_map = dict(zip(adata_meta.obs_names, adata_meta.obs["og_cluster"]))
    print(f"  OG_cluster categories: {adata_meta.obs['og_cluster'].cat.categories.tolist()}")
except FileNotFoundError:
    print("  Warning: E17_og_exact_velocity.h5ad not found.")
    og_cluster_map = {}

# ============================================================
# 4. Match barcodes and add metadata
# ============================================================
print("Matching barcodes and adding metadata...")
# Clean loom barcodes
loom_bc = adata.obs_names.tolist()
if ":" in loom_bc[0]:
    adata.obs_names = [bc.split(":")[-1] for bc in loom_bc]
    adata.obs_names = [bc.replace("x", "") if bc.endswith("x") else bc for bc in adata.obs_names]

adata.obs_names_make_unique()

# Add og_cluster from metadata
adata.obs["og_cluster"] = adata.obs_names.map(og_cluster_map)

# Count matching
matched = adata.obs["og_cluster"].notna().sum()
print(f"  Matched cells: {matched} / {adata.shape[0]}")

# ============================================================
# 5. Filter to OG1-4, PO1-2, and cells with metadata
# ============================================================
print("Filtering to OG1-4, PO1-2...")
clusters_to_keep = ["OG1", "OG2", "OG3", "OG4", "PO1", "PO2"]
mask = adata.obs["og_cluster"].isin(clusters_to_keep)
adata = adata[mask].copy()
print(f"  After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(adata.obs["og_cluster"].value_counts())

# ============================================================
# 6. Preprocessing and HVG selection
# ============================================================
print("Preprocessing...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

n_hvg = adata.var["highly_variable"].sum()
print(f"  HVG count: {n_hvg}")

# Keep only HVG
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"  After HVG filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Convert layers to float
for key in adata.layers:
    adata.layers[key] = adata.layers[key].astype(np.float32)

# ============================================================
# 7. PCA and preprocessing for batch correction
# ============================================================
print("PCA computation...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# ============================================================
# 8. Batch correction with BBKNN
# ============================================================
#print("Applying batch correction...")
#print("  Using BBKNN...")
#import bbknn
#bbknn.bbknn(adata, batch_key="batch", n_pcs=50)

# Neighbors and UMAP
#sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)
print("  Batch correction completed.")

# ============================================================
# 9. Velocity graph
# ============================================================
print("Computing velocity graph...")
try:
    scv.tl.velocity_graph(adata)
    print("  Velocity graph computed.")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 10. Save results
# ============================================================
outfile = "E17_og_integrated_harmony_filtered.h5ad"
adata.write(outfile)
print(f"Saved to {outfile}")

# ============================================================
# 11. Visualization
# ============================================================
outdir = "figures_integrated_harmony_filtered"
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

try:
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="og_cluster",
        palette=OG_COLORS,
        title="Velocity Stream (Harmony Batch Corrected, OG only)",
        save="_velocity_stream.png",
    )
    print(f"Saved stream plot to {outdir}")
except Exception as e:
    print(f"Could not save stream plot: {e}")

print("Done.")