"""
OG1-4 サブセット RNA Velocity 解析（確定版）
GSE163693_E15_17_composite_metadata.csv.gz の正確なクラスターラベルを使用
"""

import os
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

warnings.filterwarnings("ignore")
scv.settings.verbosity = 3
OUTDIR = "figures_og_exact"
sc.settings.figdir = OUTDIR
os.makedirs(OUTDIR, exist_ok=True)

WORKDIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORKDIR)

# ============================================================
# 1. メタデータ読み込み → E17 OG1-4 バーコードを取得
# ============================================================
print("Loading metadata ...")
meta = pd.read_csv("GSE163693_E15_17_composite_metadata.csv.gz", index_col=0)

meta_e17 = meta[meta["orig.ident"] == "E17hiseq"].copy()
print(f"  E17 cells in metadata: {len(meta_e17)}")
print(f"  Cluster distribution:\n{meta_e17['Cluster'].value_counts().to_string()}")

# OG1-4 に絞る
og_labels = ["OG1", "OG2", "OG3", "OG4"]
meta_og = meta_e17[meta_e17["Cluster"].isin(og_labels)].copy()
print(f"\n  OG1-4 cells in E17 metadata: {len(meta_og)}")
print(meta_og["Cluster"].value_counts().to_string())

# バーコードの末尾 _1 を除去して velocity AnnData 形式に合わせる
# 例: AAACCTGAGACAATAC-1_1 → AAACCTGAGACAATAC-1
meta_og.index = meta_og.index.str.replace(r"_\d+$", "", regex=True)
og_barcodes = set(meta_og.index)
print(f"\n  OG barcodes after stripping suffix: {len(og_barcodes)} (e.g. {list(og_barcodes)[:3]})")

# ============================================================
# 2. velocity AnnData を読み込み、OG1-4 にサブセット
# ============================================================
print("\nLoading E17_suture_velocity.h5ad ...")
adata_full = sc.read_h5ad("E17_suture_velocity.h5ad")
print(f"  Full AnnData: {adata_full.shape[0]} cells x {adata_full.shape[1]} genes")

common = list(set(adata_full.obs_names) & og_barcodes)
print(f"  Matching OG barcodes: {len(common)}")

adata_og = adata_full[common].copy()

# メタデータのクラスターラベルを付与
adata_og.obs["og_cluster"] = meta_og.loc[adata_og.obs_names, "Cluster"]
print(f"\nOG subset: {adata_og.shape[0]} cells")
print(adata_og.obs["og_cluster"].value_counts().to_string())

# ============================================================
# 3. 全体 UMAP 上で OG1-4 の位置を確認
# ============================================================
adata_full.obs["og_cluster"] = "other"
adata_full.obs.loc[common, "og_cluster"] = meta_og.loc[common, "Cluster"]

sc.pl.umap(
    adata_full,
    color="og_cluster",
    palette=["#e41a1c", "#ff7f00", "#4daf4a", "#984ea3", "#cccccc"],
    title="OG1-4 location in full UMAP",
    save="_og_in_full_umap.png",
)
print(f"  Saved {OUTDIR}/umap_og_in_full_umap.png")

# ============================================================
# 4. loom ファイルから spliced/unspliced を読み込んでマージ
# ============================================================
def read_loom_as_anndata(filename):
    with loompy.connect(filename, "r") as ds:
        X = sp.csc_matrix(ds[:, :]).T
        layers = {}
        for key in ds.layers.keys():
            if key == "":
                continue
            layers[key] = sp.csc_matrix(ds.layers[key][:, :]).T
        obs_names = pd.Index(ds.col_attrs["CellID"][:])
        var_names = pd.Index(ds.row_attrs["Gene"][:])
        ldata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=var_names),
        )
        for key, mat in layers.items():
            ldata.layers[key] = mat
    return ldata

import glob
loom_files = sorted(glob.glob("loom_output/*.loom"))
loom_list = []
for f in loom_files:
    print(f"  Loading loom: {f}")
    ldata = read_loom_as_anndata(f)
    ldata.var_names_make_unique()
    loom_list.append(ldata)

ldata_merged = ad.concat(loom_list, join="outer") if len(loom_list) > 1 else loom_list[0]
loom_barcodes = ldata_merged.obs_names.tolist()
if ":" in loom_barcodes[0]:
    ldata_merged.obs_names = [bc.split(":")[-1] for bc in ldata_merged.obs_names]
    ldata_merged.obs_names = [
        bc[:-1] if bc.endswith("x") else bc for bc in ldata_merged.obs_names
    ]
ldata_merged.obs_names_make_unique()
ldata_merged.var_names_make_unique()

adata_og_vel = scv.utils.merge(adata_og, ldata_merged)
print(f"  OG merged with loom: {adata_og_vel.shape[0]} cells x {adata_og_vel.shape[1]} genes")

# og_cluster ラベルを再付与（merge で失われる場合に備えて）
adata_og_vel.obs["og_cluster"] = meta_og.reindex(adata_og_vel.obs_names)["Cluster"]

# ============================================================
# 5. 前処理：OG サブセット専用に近傍グラフ・UMAP を再計算
# ============================================================
print("\nPreprocessing OG subset ...")
scv.pp.filter_and_normalize(adata_og_vel, min_shared_counts=20)
sc.pp.highly_variable_genes(adata_og_vel, n_top_genes=2000, flavor="seurat")

# 近傍グラフ・UMAP を OG サブセットで新規計算
sc.pp.pca(adata_og_vel)
sc.pp.neighbors(adata_og_vel, n_pcs=30, n_neighbors=30)
sc.tl.umap(adata_og_vel)

# scVelo moments
scv.pp.moments(adata_og_vel, n_pcs=30, n_neighbors=30)

# og_cluster をカテゴリカル型に変換して色指定
og_palette = ["#e41a1c", "#ff7f00", "#4daf4a", "#984ea3"]  # OG1,2,3,4
adata_og_vel.obs["og_cluster"] = pd.Categorical(
    adata_og_vel.obs["og_cluster"],
    categories=["OG1", "OG2", "OG3", "OG4"],
)
adata_og_vel.uns["og_cluster_colors"] = og_palette

# UMAP 確認（OG クラスターラベル）
sc.pl.umap(
    adata_og_vel,
    color="og_cluster",
    palette=["#e41a1c", "#ff7f00", "#4daf4a", "#984ea3"],
    title="OG1-4 subset UMAP",
    save="_og_exact_clusters.png",
)
print(f"  Saved {OUTDIR}/umap_og_exact_clusters.png")

# spliced/unspliced の比率確認
scv.pl.proportions(adata_og_vel, groupby="og_cluster",
                   save="_og_proportions.png")

# ============================================================
# 6. Dynamical velocity
# ============================================================
print("\nRecovering dynamics (OG1-4, dynamical model) ...")
scv.tl.recover_dynamics(adata_og_vel, n_jobs=4)
scv.tl.velocity(adata_og_vel, mode="dynamical")
scv.tl.velocity_graph(adata_og_vel)

# velocity ストリーム（OG クラスター色分け）
scv.pl.velocity_embedding_stream(
    adata_og_vel,
    basis="umap",
    color="og_cluster",
    title="OG1-4 RNA velocity stream",
    save="_og_exact_stream.png",
)
scv.pl.velocity_embedding(
    adata_og_vel,
    basis="umap",
    arrow_length=3,
    arrow_size=2,
    color="og_cluster",
    save="_og_exact_arrows.png",
)
print(f"  Saved velocity stream and arrows")

# ============================================================
# 7. Latent time（OG1 を根として）
# ============================================================
scv.tl.latent_time(adata_og_vel)

scv.pl.scatter(
    adata_og_vel,
    basis="umap",
    color="latent_time",
    color_map="gnuplot",
    title="OG1-4 Latent time",
    save="_og_exact_latent_time.png",
)

# OG クラスターごとの latent time 分布
fig, ax = plt.subplots(figsize=(7, 4))
for og, col in [("OG1", "#e41a1c"), ("OG2", "#ff7f00"),
                ("OG3", "#4daf4a"), ("OG4", "#984ea3")]:
    cells = adata_og_vel.obs[adata_og_vel.obs["og_cluster"] == og]["latent_time"]
    if len(cells) > 0:
        cells.hist(bins=30, alpha=0.6, label=og, color=col, ax=ax, density=True)
ax.set_xlabel("Latent time")
ax.set_ylabel("Density")
ax.set_title("Latent time distribution per OG cluster")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTDIR}/og_latent_time_hist.png", dpi=150)
plt.close()
print(f"  Saved {OUTDIR}/og_latent_time_hist.png")

# ============================================================
# 8. OG マーカー遺伝子の latent time に沿った発現
# ============================================================
key_genes_priority = [
    "Erg", "Pthlh", "Six2",        # OG1
    "Lef1", "Inhba",               # OG2
    "Mmp13", "Podnl1",             # OG3
    "Ifitm5", "Dmp1", "Sost",      # OG4
    "Runx2", "Sp7",                # 骨芽細胞全般
]
available_keys = [g for g in key_genes_priority if g in adata_og_vel.var_names]
print(f"  Available key genes for plots: {available_keys}")

if available_keys:
    scv.pl.scatter(
        adata_og_vel,
        x="latent_time",
        y=available_keys,
        color="og_cluster",
        ncols=4,
        frameon=False,
        save="_og_exact_markers_vs_time.png",
    )
    sc.pl.umap(
        adata_og_vel,
        color=available_keys,
        ncols=4,
        save="_og_exact_key_genes.png",
    )

# ============================================================
# 9. ドライバー遺伝子
# ============================================================
scv.tl.rank_velocity_genes(adata_og_vel, groupby="og_cluster", min_corr=0.3)
df_driver = pd.DataFrame(adata_og_vel.uns["rank_velocity_genes"]["names"])
df_driver.to_csv("og_exact_velocity_driver_genes.csv", index=False)
print("\nDriver genes (top 10 per OG cluster):")
print(df_driver.head(10).to_string())

# ============================================================
# 10. 保存
# ============================================================
adata_og_vel.write("E17_og_exact_velocity.h5ad")

print("\n=== DONE ===")
print(f"  OG cells: {adata_og_vel.n_obs}")
print(f"    {adata_og_vel.obs['og_cluster'].value_counts().to_string()}")
print(f"  Figures: {OUTDIR}/")
print(f"  AnnData: E17_og_exact_velocity.h5ad")
print(f"  Driver genes: og_exact_velocity_driver_genes.csv")
