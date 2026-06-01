"""
HVGのみを使ってHarmonyで統合するスクリプト（メモリ効率版）

- loom_output/ のLoomファイルを読み込み
- HVG（高変動遺伝子）のみに削減
- Harmonyでバッチ補正
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

try:
    from harmony import harmonize
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False
    raise ImportError("harmonyphpy is required. Install with: pip install harmonypy")

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

def read_loom_as_anndata(filename):
    """Read a velocyto loom file into AnnData."""
    import loompy
    with loompy.connect(filename, "r") as ds:
        # spliced layer をX に設定
        if "spliced" in ds.layers:
            X = ds.layers["spliced"][:, :].T  # cells x genes
        else:
            X = np.zeros((ds.shape[1], ds.shape[0]))
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

    # バッチ情報をファイル名から推定
    basename = os.path.basename(f)
    batch_id = basename.split("batch_")[1].split("_")[0] if "batch_" in basename else basename.split(".")[0]
    ldata.obs["batch"] = f"batch_{batch_id}"

    loom_list.append(ldata)
    print(f"  {ldata.n_obs} cells x {ldata.n_vars} genes")

# ============================================================
# 2. 統合（concat）
# ============================================================
print("\nIntegrating loom files...")
if len(loom_list) == 1:
    adata = loom_list[0]
else:
    adata = ad.concat(loom_list, join="outer")

print(f"Integrated data: {adata.n_obs} cells x {adata.n_vars} genes")

# ============================================================
# 3. メモリ効率的な前処理
# ============================================================
print("\nPreprocessing...")

# Layer型を修正
for key in adata.layers:
    adata.layers[key] = adata.layers[key].astype(np.float32)

# QC フィルタリング
print("Filtering cells and genes...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

# 正規化とHVG検出
print("Normalizing and detecting HVGs...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)

# HVGのみに削減
adata = adata[:, adata.var.highly_variable].copy()
print(f"After HVG selection: {adata.n_obs} cells x {adata.n_vars} genes")

# ============================================================
# 4. PCA計算
# ============================================================
print("Computing PCA...")
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50)

# ============================================================
# 5. Harmony でバッチ補正
# ============================================================
print("\nApplying Harmony batch correction...")
adata = harmonize(adata, key="batch", max_iter_harmony=10)
print("Harmony completed.")

# ============================================================
# 6. UMAP 計算
# ============================================================
print("Computing neighbors and UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

# ============================================================
# 7. 保存
# ============================================================
output_file = "integrated_loom_harmony.h5ad"
print(f"\nSaving to {output_file}...")
adata.write(output_file)
print(f"Saved integrated data to {output_file}")

print("Done.")
