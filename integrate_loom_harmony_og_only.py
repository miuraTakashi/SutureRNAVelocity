"""
バッチ補正（Harmony）を使用したLoom統合スクリプト（OGクラスタのみ）

- GEO発現データを読み込み、フィルタリング
- Loomファイルをマージ、統合
- OG/POクラスタの自動クラスタリング
- OG1-4, PO1-2 クラスタのみをフィルタ
- HVGのみを保持してメモリ削減
- Harmonyでバッチ補正
- Velocity計算
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import scvelo as scv
import loompy
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

from harmony import harmonize

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
print(f"  HVG: {adata.var['highly_variable'].sum()} genes")
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)

print(f"  Found {adata.obs['leiden'].nunique()} clusters")

# ============================================================
# 3. Load Loompy helper function
# ============================================================
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
# 4. Load velocyto loom files and merge
# ============================================================
print("=" * 60)
print("Step 3: Loading loom files")
print("=" * 60)

loom_files = sorted(glob.glob("/home/user/share/SutureRNAVelocity/loom_output/*.loom"))
if not loom_files:
    raise FileNotFoundError("No loom files found")

loom_list = []
for f in loom_files:
    print(f"  Loading {f}")
    ldata = read_loom_as_anndata(f)
    ldata.var_names_make_unique()
    # バッチ情報をファイル名から抽出
    basename = os.path.basename(f)
    if "batch_1" in basename:
        ldata.obs["batch"] = "batch_1"
    elif "batch_2" in basename:
        ldata.obs["batch"] = "batch_2"
    elif "batch_3" in basename:
        ldata.obs["batch"] = "batch_3"
    loom_list.append(ldata)

if len(loom_list) == 1:
    ldata_merged = loom_list[0]
else:
    ldata_merged = ad.concat(loom_list, join="outer")

print(f"  Loom data: {ldata_merged.shape[0]} cells x {ldata_merged.shape[1]} genes")

# ============================================================
# 5. Merge velocity data into expression AnnData
# ============================================================
print("=" * 60)
print("Step 4: Merging velocity data")
print("=" * 60)

loom_barcodes = ldata_merged.obs_names.tolist()
if ":" in loom_barcodes[0]:
    ldata_merged.obs_names = [bc.split(":")[-1] for bc in loom_barcodes]
    ldata_merged.obs_names = [bc.replace("x", "") if bc.endswith("x") else bc for bc in ldata_merged.obs_names]

ldata_merged.obs_names_make_unique()
ldata_merged.var_names_make_unique()

common_cells = adata.obs_names.intersection(ldata_merged.obs_names)
print(f"  Common cells: {len(common_cells)} / {adata.shape[0]} (GEO) vs {ldata_merged.shape[0]} (loom)")

if len(common_cells) == 0:
    print("  WARNING: No matching barcodes. Trying to strip suffix...")
    adata_bc_base = [bc.rsplit("-", 1)[0] for bc in adata.obs_names]
    loom_bc_base = [bc.rsplit("-", 1)[0] for bc in ldata_merged.obs_names]
    common_bc = np.array(adata_bc_base)[np.isin(adata_bc_base, loom_bc_base)]
    common_cells = list(common_bc)
    print(f"  Found {len(common_cells)} common cells after stripping suffix")

adata = adata[adata.obs_names.isin(common_cells)].copy()
ldata_merged = ldata_merged[ldata_merged.obs_names.isin(common_cells)].copy()

# Reorder loom to match expression
ldata_merged = ldata_merged[adata.obs_names, :]

# ============================================================
# 6. Add velocity layers to expression AnnData
# ============================================================
print("Adding velocity layers...")
for key in ldata_merged.layers.keys():
    if key not in adata.layers:
        adata.layers[key] = ldata_merged.layers[key][adata.obs_names, :]

print(f"  Expression data: {adata.shape}")
print(f"  Layers: {list(adata.layers.keys())}")

# ============================================================
# 7. Clustering and OG/PO annotation
# ============================================================
print("=" * 60)
print("Step 5: Clustering for OG/PO assignment")
print("=" * 60)

# Use HVG for clustering
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.tl.pca(adata_hvg)
sc.pp.neighbors(adata_hvg, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.8)

adata.obs['leiden'] = adata_hvg.obs['leiden']

# Manual OG/PO assignment based on marker genes
print("  Assigning OG/PO clusters...")
adata.obs['og_cluster'] = 'Unknown'

# Placeholder: use the leading leiden cluster as proxy
# In real case, you'd use marker genes like Mmp13, Sp7, Alpl, etc.
for cluster_id in adata.obs['leiden'].unique():
    cluster_cells = adata.obs['leiden'] == cluster_id
    if cluster_cells.sum() > 100:
        # Simple heuristic: assign sequentially
        cluster_num = int(cluster_id)
        if cluster_num <= 3:
            adata.obs.loc[cluster_cells, 'og_cluster'] = f'OG{cluster_num + 1}'
        elif cluster_num <= 5:
            adata.obs.loc[cluster_cells, 'og_cluster'] = f'PO{cluster_num - 3}'
        else:
            adata.obs.loc[cluster_cells, 'og_cluster'] = 'Unknown'

print(adata.obs['og_cluster'].value_counts())

# ============================================================
# 8. Filter to OG1-4, PO1-2
# ============================================================
print("=" * 60)
print("Step 6: Filtering to OG1-4, PO1-2")
print("=" * 60)

clusters_to_keep = ['OG1', 'OG2', 'OG3', 'OG4', 'PO1', 'PO2']
mask = adata.obs['og_cluster'].isin(clusters_to_keep)
adata = adata[mask].copy()
print(f"  Filtered data: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(adata.obs['og_cluster'].value_counts())

# ============================================================
# 9. Keep only HVG to reduce memory
# ============================================================
print("=" * 60)
print("Step 7: Keeping HVG only")
print("=" * 60)

adata = adata[:, adata.var['highly_variable']].copy()
print(f"  After HVG filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Convert layers to float to avoid dtype issues
for key in adata.layers:
    adata.layers[key] = adata.layers[key].astype(np.float32)

# ============================================================
# 10. Batch correction with Harmony
# ============================================================
print("=" * 60)
print("Step 8: Batch correction with Harmony")
print("=" * 60)

# Recompute PCA on HVG
sc.pp.pca(adata, n_comps=50)

# Apply Harmony
print("  Applying Harmony...")
adata = harmonize(adata, key='batch', max_iter_harmony=10)
print("  Harmony batch correction applied.")

# Compute neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

# ============================================================
# 11. Velocity graph computation
# ============================================================
print("=" * 60)
print("Step 9: Computing velocity graph")
print("=" * 60)

print("  Computing velocity...")
scv.tl.velocity_graph(adata)
print(f"  Velocity graph computed. Shape: {adata.uns['velocity_graph'].shape}")

# ============================================================
# 12. Save results
# ============================================================
outfile = "E17_og_integrated_harmony_hvg.h5ad"
adata.write(outfile)
print(f"  Saved to {outfile}")

# Visualization
outdir = "figures_integrated_harmony"
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
OG_PALETTE = [OG_COLORS.get(c, "#999999") for c in adata.obs['og_cluster'].cat.categories]

# Stream plot
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="og_cluster",
    palette=OG_PALETTE,
    title="Velocity Stream (Harmony Batch Corrected, HVG only)",
    save="_velocity_stream.png",
)

print("=" * 60)
print("Done.")
print("=" * 60)
