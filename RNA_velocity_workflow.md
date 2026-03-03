# RNA Velocity 解析ワークフロー

## 概要

GEO: GSE163693 のマウス頭蓋縫合 (cranial suture) scRNA-seq データ（E17, 3バッチ）を用いて RNA velocity 解析を行う。
10x Chromium の BAM ファイルから spliced/unspliced カウント行列を生成し、scVelo で RNA velocity を推定する。

## データ構成

```
SutureRNAVelocity/
├── ena_runs.tsv                        # ENA ダウンロード情報
├── gencode.vM25.annotation.gtf         # マウスゲノムアノテーション (GRCm38)
├── GSE163693_RAW/                      # GEO 発現行列（Cell Ranger 出力）
│   ├── GSM4983998_E15_*                #   E15 サンプル
│   └── GSM4983999_E17_*                #   E17 サンプル
└── downloads/                          # ENA からダウンロードした BAM
    ├── SRR13288596/E17_batch_1_possorted_genome_bam.bam{,.bai}
    ├── SRR13288597/E17_batch_2_possorted_genome_bam.bam{,.bai}
    └── SRR13288598/E17_batch_3_possorted_genome_bam.bam{,.bai}
```

---

## Step 0: 環境構築

### conda 環境の作成

```bash
conda create -n velocity python=3.10 -y
conda activate velocity
```

### パッケージのインストール

```bash
pip install velocyto scvelo scanpy loompy anndata matplotlib
# velocyto CLI には samtools も必要
conda install -c bioconda samtools -y
```

---

## Step 1: リピートマスク GTF の準備

velocyto はリピート領域をマスクするための GTF ファイルを推奨している。
UCSC Table Browser から mm10 の RepeatMasker アノテーションを取得する。

```bash
# UCSC から mm10 リピートマスク GTF をダウンロード
curl -L -o mm10_rmsk.gtf.gz \
  "https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=&clade=mammal&org=Mouse&db=mm10&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf.gz"
gunzip mm10_rmsk.gtf.gz
```

> **注意**: UCSC Table Browser からの自動取得がうまくいかない場合は、ブラウザで手動ダウンロードするか、以下の代替方法を使う。

```bash
# 代替: UCSC の rmsk テーブルから直接取得
wget "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"
gunzip rmsk.txt.gz
# rmsk.txt → GTF 形式に変換
awk 'BEGIN{OFS="\t"} {print $6,"rmsk",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$11"\";"}' rmsk.txt > mm10_rmsk.gtf
```

---

## Step 2: セルバーコードの準備

10x Chromium v2/v3 のホワイトリストバーコードを使用する。
バーコードは Cell Ranger のインストールに含まれているが、単体でも取得可能。

```bash
# 10x Chromium v3 バーコード（Cell Ranger に同梱）
# Cell Ranger がインストール済みの場合:
#   cellranger-x.x.x/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

# または GEO に含まれるバーコードファイルを使用
zcat GSE163693_RAW/GSM4983999_E17_barcodes.tsv.gz > E17_barcodes.tsv
```

> **ポイント**: GEO の `barcodes.tsv.gz` には実際にフィルタリングされたバーコードのみが含まれているため、これを使用すれば該当サンプルの有効セルに限定した解析が可能になる。ただし、3バッチ統合 BAM の場合はバーコードのバッチ情報（サフィックス `-1`, `-2` 等）に注意が必要。

---

## Step 3: velocyto run10x による Loom ファイル生成

各バッチの BAM ファイルから spliced/unspliced/ambiguous カウント行列を含む `.loom` ファイルを生成する。

```bash
ANNOTATION="gencode.vM25.annotation.gtf"
RMSK="mm10_rmsk.gtf"   # Step 1 で作成

# Batch 1
velocyto run \
  -b E17_barcodes.tsv \
  -o loom_output \
  -m $RMSK \
  downloads/SRR13288596/E17_batch_1_possorted_genome_bam.bam \
  $ANNOTATION

# Batch 2
velocyto run \
  -b E17_barcodes.tsv \
  -o loom_output \
  -m $RMSK \
  downloads/SRR13288597/E17_batch_2_possorted_genome_bam.bam \
  $ANNOTATION

# Batch 3
velocyto run \
  -b E17_barcodes.tsv \
  -o loom_output \
  -m $RMSK \
  downloads/SRR13288598/E17_batch_3_possorted_genome_bam.bam \
  $ANNOTATION
```

> **所要時間の目安**: BAM ファイルのサイズ（20〜39 GB）に依存するが、1バッチあたり数時間〜半日程度。
> **メモリ**: 16 GB 以上推奨。

---

## Step 4: Loom ファイルの統合と前処理 (Python)

```python
import scvelo as scv
import scanpy as sc
import loompy
import numpy as np

scv.settings.verbosity = 3
scv.settings.presenter_view = True

# --- 4.1 Loom ファイルの読み込み ---
loom_files = [
    "loom_output/E17_batch_1_possorted_genome_bam.loom",
    "loom_output/E17_batch_2_possorted_genome_bam.loom",
    "loom_output/E17_batch_3_possorted_genome_bam.loom",
]

# 複数 loom を結合
ldata = scv.read(loom_files[0])
for f in loom_files[1:]:
    tmp = scv.read(f)
    ldata = ldata.concatenate(tmp)

# --- 4.2 GEO の発現行列との統合（任意） ---
# 既存のクラスタ情報や UMAP 座標を活用する場合
adata = sc.read_10x_mtx(
    "GSE163693_RAW/",
    prefix="GSM4983999_E17_",
    var_names="gene_ids",
)
# バーコードの対応付け
# velocyto 出力のバーコードは "sample:barcode" 形式の場合がある
# 必要に応じてバーコードを整形してマージ
scv.utils.merge(adata, ldata)
```

---

## Step 5: scVelo による RNA Velocity 解析

```python
# --- 5.1 フィルタリングと正規化 ---
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# --- 5.2 Velocity の推定 ---
# Stochastic モデル（標準）
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# --- 5.3 可視化 ---
scv.pl.velocity_embedding_stream(adata, basis="umap", save="velocity_stream.png")
scv.pl.velocity_embedding(adata, basis="umap", arrow_length=3, arrow_size=2, save="velocity_arrows.png")
```

---

## Step 6: Dynamical モデルによる高精度解析（推奨）

```python
# Dynamical モデルは遺伝子ごとに転写・スプライシング・分解の速度定数を推定する
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Latent time の推定
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", save="latent_time.png")

# ドライバー遺伝子のランキング
scv.tl.rank_velocity_genes(adata, groupby="clusters", min_corr=0.3)
df = scv.DataFrame(adata.uns["rank_velocity_genes"]["names"])
df.head(10)
```

---

## Step 7: 結果の保存

```python
adata.write("E17_suture_velocity.h5ad")

# 主要な図の出力
scv.pl.velocity_embedding_stream(adata, basis="umap", save="final_velocity_stream.pdf")
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", save="final_latent_time.pdf")
```

---

## トラブルシューティング

| 問題 | 対処法 |
|------|--------|
| velocyto が遅い / メモリ不足 | `samtools sort -t CB` で BAM をバーコード順にソートしてから実行 |
| Loom ファイルのバーコードが一致しない | バーコードのプレフィックス/サフィックスを確認し `adata.obs_names` を整形 |
| Velocity が計算されない遺伝子が多い | `min_shared_counts` を下げる、`n_top_genes` を増やす |
| UMAP 上で矢印が無秩序 | 前処理パラメータ調整、dynamical モデルの使用を検討 |
| `No spliced/unspliced layers` エラー | Loom ファイルの統合時に layer が正しくマージされているか確認 |

---

## 参考リンク

- [scVelo documentation](https://scvelo.readthedocs.io/)
- [velocyto documentation](http://velocyto.org/)
- [La Manno et al., 2018 (RNA velocity)](https://doi.org/10.1038/s41586-018-0414-6)
- [Bergen et al., 2020 (scVelo)](https://doi.org/10.1038/s41587-020-0591-3)
- [GSE163693 (GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163693)
