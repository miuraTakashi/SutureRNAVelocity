"""
Loomファイルをバッチ補正を行った上で統合するスクリプト

- loom_output/ のLoomファイルを読み込み
- バッチ補正（HarmonyまたはBBKNN）を適用
- 統合してh5adファイルを保存
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore")

# Batch correction 用
try:
    from harmonypy import run_harmony
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False
    print("Warning: harmonypy not installed. Will use BBKNN.")

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

def read_loom_as_anndata(filename):
    """Read a velocyto loom file into AnnData."""
    import loompy
    with loompy.connect(filename, "r") as ds:
        # X を spliced layer に設定
        if "spliced" in ds.layers:
            X = ds.layers["spliced"][:, :].T  # cells x genes
        else:
            X = np.zeros((ds.shape[1], ds.shape[0]))  # cells x genes
        layers = {}
        for key in ds.layers.keys():
            if key == "":
                continue
            layers[key] = ds.layers[key][:, :].T
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

# ============================================================
# 1. Loomファイルの読み込み
# ============================================================
loom_dir = "/home/user/share/SutureRNAVelocity/loom_output"
loom_files = sorted(glob.glob(os.path.join(loom_dir, "*.loom")))
if not loom_files:
    raise FileNotFoundError(f"No loom files found in {loom_dir}/")

print(f"Found {len(loom_files)} loom files:")
for f in loom_files:
    print(f"  {f}")

loom_list = []
for f in loom_files:
    print(f"Loading {f}...")
    ldata = read_loom_as_anndata(f)
    ldata.var_names_make_unique()

    # バッチ情報をファイル名から推定（例: SRR ID）
    # ファイル名にSRRが含まれていると仮定
    basename = os.path.basename(f)
    if "SRR" in basename:
        batch_id = basename.split("SRR")[1].split("_")[0]  # 例: SRR123456
        ldata.obs["batch"] = f"SRR{batch_id}"
    else:
        ldata.obs["batch"] = basename.split(".")[0]  # ファイル名をバッチとして

    loom_list.append(ldata)

# ============================================================
# 2. 統合（concat）
# ============================================================
print("Integrating loom files...")
if len(loom_list) == 1:
    adata_integrated = loom_list[0]
else:
    adata_integrated = ad.concat(loom_list, join="outer")

print(f"Integrated data: {adata_integrated.shape[0]} cells x {adata_integrated.shape[1]} genes")

# 前処理（PCA用）
print("Preprocessing for PCA...")
print("  Filtering cells and genes...")
sc.pp.filter_cells(adata_integrated, min_genes=200)
sc.pp.filter_genes(adata_integrated, min_cells=3)
print(f"  After filtering: {adata_integrated.n_obs} cells x {adata_integrated.n_vars} genes")

print("  Normalizing and log-transforming...")
sc.pp.normalize_total(adata_integrated, target_sum=1e4)
sc.pp.log1p(adata_integrated)

print("  Computing highly variable genes...")
sc.pp.highly_variable_genes(adata_integrated, min_mean=0.0125, max_mean=3, min_disp=0.5)

print("  Scaling...")
sc.pp.scale(adata_integrated, max_value=10)
print("Preprocessing complete.")

# ============================================================
# 3. バッチ補正
# ============================================================
print("\n=== Batch Correction ===")
if "batch" in adata_integrated.obs and adata_integrated.obs["batch"].nunique() > 1:
    print(f"Applying batch correction (batches: {adata_integrated.obs['batch'].nunique()})...")
    batch_key = "batch"

    if HAS_HARMONY:
        print("Using Harmony for batch correction...")
        try:
            print("  Computing PCA...")
            sc.pp.pca(adata_integrated, n_comps=30, svd_solver='arpack')
            print("  Applying Harmony...")
            from harmony import harmonize
            adata_corrected = harmonize(
                adata_integrated,
                key=batch_key,
                max_iter_harmony=10,
                verbose=False
            )
            print("  Harmony batch correction completed.")
        except Exception as e:
            print(f"  Harmony error: {e}")
            print("  Falling back to BBKNN...")
            import bbknn
            sc.pp.pca(adata_integrated, n_comps=30, svd_solver='arpack')
            bbknn.bbknn(adata_integrated, batch_key=batch_key, n_pcs=30)
            adata_corrected = adata_integrated.copy()
    else:
        print("Using BBKNN for batch correction...")
        print("  Computing PCA...")
        sc.pp.pca(adata_integrated, n_comps=30, svd_solver='arpack')
        print("  Applying BBKNN...")
        import bbknn
        bbknn.bbknn(adata_integrated, batch_key=batch_key, n_pcs=30)
        adata_corrected = adata_integrated.copy()

    # PCA/UMAP再計算
    print("  Computing neighbors and UMAP...")
    sc.pp.neighbors(adata_corrected, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_corrected)

    adata_integrated = adata_corrected
    print("Batch correction applied.")
else:
    print("No batch correction needed (single batch or no batch info).")

# ============================================================
# 4. 保存
# ============================================================
output_file = "integrated_velocity_batch_corrected.h5ad"
adata_integrated.write(output_file)
print(f"Saved integrated data to {output_file}")

print("Done.")