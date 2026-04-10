# SutureRNAVelocity

マウス頭蓋縫合 scRNA-seq データ（E17.5）の RNA velocity 解析パイプライン。

- **データソース**: GEO [GSE163693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163693)
- **論文**: Farmer et al. 2021, *Nature Communications* 12, 4797 — [DOI](https://doi.org/10.1038/s41467-021-24917-9)
- **サンプル**: E17.5 マウス冠状縫合、3 バッチ（SRR13288596–SRR13288598）
- **ゲノムアノテーション**: Gencode vM25（GRCm38/mm10）

---

## 環境構築

```bash
# Python 仮想環境の作成
python3 -m venv .venv
source .venv/bin/activate

# 基本パッケージ
pip install numpy scipy cython pysam setuptools wheel
pip install --no-build-isolation velocyto

# 解析パッケージ
pip install scanpy scvelo anndata loompy leidenalg igraph pandas matplotlib
pip install cellrank
```

---

## 解析ワークフロー

### Step 1: BAM ファイルのダウンロード

```bash
aria2c --input-file=aria2c_download.txt \
       --max-concurrent-downloads=3 \
       --split=8 \
       --max-connection-per-server=8 \
       --continue=true
```

- ダウンロード先: `downloads/SRR13288596/`, `SRR13288597/`, `SRR13288598/`
- BAM ファイル + BAI インデックスファイルを取得

### Step 2: GTF アノテーションの準備

```bash
# chr プレフィックスを除去（velocyto の要件）
sed 's/^chr//' gencode.vM25.annotation.gtf > gencode.vM25.nochr.raw.gtf

# transcript_id を含む行のみ抽出
awk -F '\t' 'BEGIN{OFS="\t"} NF>=9 && $9 ~ /transcript_id/ \
    {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
    gencode.vM25.nochr.raw.gtf > gencode.vM25.nochr.gtf
```

### Step 3: velocyto による spliced/unspliced カウント

```bash
bash run_velocyto.sh
```

- 出力: `loom_output/E17_batch_1_*.loom`, `batch_2`, `batch_3`
- バーコードファイル: `E17_barcodes.tsv`
- velocyto 内部バグへのパッチ適用済み（`counter.py` の `nth_line` 初期化）

### Step 4: 全細胞 RNA velocity 解析（scVelo）

```bash
source .venv/bin/activate
python run_scvelo.py
```

**処理内容:**
1. GEO 発現行列（`GSE163693_RAW/`）読み込み
2. loom ファイルと結合（spliced/unspliced マージ）
3. `scv.pp.filter_and_normalize` → highly variable genes 選択
4. `scv.tl.recover_dynamics` → dynamical モデルで velocity 推定
5. Leiden クラスタリング（19 クラスター）+ UMAP
6. Latent time・ドライバー遺伝子の算出

**主な出力:**
| ファイル | 内容 |
|---|---|
| `E17_suture_velocity.h5ad` | 全細胞 velocity 解析済みデータ |
| `velocity_driver_genes.csv` | クラスター別ドライバー遺伝子 |
| `figures/umap_clusters.png` | 全細胞 UMAP（Leiden クラスター） |
| `figures/scvelo_velocity_stream_dynamical.png` | velocity ストリーム |
| `figures/scvelo_latent_time.png` | Latent time |

---

### Step 5: OG1–OG4 サブセット解析

論文（Farmer et al. 2021）で定義された骨形成細胞群（OG1–OG4）のみを抽出して velocity を再計算。

#### クラスター定義（Farmer et al. 2021）

| クラスター | 細胞タイプ | マーカー遺伝子 |
|:-:|---|---|
| **OG1** | 骨形成前駆細胞（最初期） | *Erg*, *Pthlh*, *Six2* |
| **OG2** | 縫合部プレ骨芽細胞 | *Lef1*, *Inhba* |
| **OG3** | 骨膜プレ骨芽細胞 | *Mmp13*, *Podnl1* |
| **OG4** | 成熟骨芽細胞 | *Ifitm5*, *Dmp1*, *Sost* |

#### バーコード同定

`GSE163693_E15_17_composite_metadata.csv.gz` の `Cluster` 列を使用。E17.5 の OG1–OG4 バーコード（計 1,772 細胞）を正確に取得。

```bash
python run_og_exact.py
```

**細胞数:**

| クラスター | 細胞数 |
|:-:|:-:|
| OG3 | 603 |
| OG1 | 548 |
| OG2 | 315 |
| OG4 | 306 |
| **合計** | **1,772** |

**主な出力:**
| ファイル | 内容 |
|---|---|
| `E17_og_exact_velocity.h5ad` | OG サブセット velocity 解析済みデータ |
| `og_exact_velocity_driver_genes.csv` | OG 別ドライバー遺伝子 |
| `figures_og_exact/velocity_stream.png` | OG サブセット velocity ストリーム |
| `figures_og_exact/latent_time.png` | Latent time |
| `figures_og_exact/proportions.png` | Spliced/unspliced 比率 |

---

### Step 6: 分化・脱分化トポロジーの解析（CellRank）

velocity ストリームがループ状に見えたため、分化方向を固定せず双方向性（脱分化 OG4→OG3 を含む）かつ、**必ず終末状態に吸収されるという仮定を置かない**形でトポロジーを評価する。

```bash
python run_og_cellrank.py
```

**PAGA 遷移強度（エッジ重み）:**

```
      OG1    OG2    OG3    OG4
OG1    —     0.772  0.752  0.059
OG2   0.772   —     0.673  0.368
OG3   0.752  0.673   —     0.125
OG4   0.059  0.368  0.125   —
```

CellRank では `VelocityKernel` + `ConnectivityKernel` により**非対称（非可逆）な遷移行列**を構築し、GPCCA によって macrostates を同定する。このとき terminal state / initial state や absorption probability は計算せず、**OG1–OG4 間の遷移構造そのもの（ループを含む）**を解釈対象とする。

さらに、`combined_kernel.transition_matrix` から OG クラスター間の正味フラックス（A→B - B→A）を計算し、UMAP 空間上のクラスター重心間に矢印として描画することで、**粗視化された分化方向ベクトル場**を得る。

| ファイル | 内容 |
|---|---|
| `figures_og_exact/og_paga_graph.png` | PAGA トポロジーグラフ |
| `figures_og_exact/cr_macrostates.png` | GPCCA マクロ状態 |
| `figures_og_exact/cr_diff_vector_field.png` | OG クラスター間の分化ベクトル場（正味フラックス） |

**主な出力:**
| ファイル | 内容 |
|---|---|
| `figures_og_exact/og_paga_graph.png` | PAGA トポロジーグラフ |
| `figures_og_exact/cr_macrostates.png` | GPCCA マクロ状態 |
| `figures_og_exact/cr_initial_terminal.png` | 初期・終末状態 |
| `figures_og_exact/cr_fate_probabilities.png` | 各細胞の運命確率 |
| `figures_og_exact/cr_fate_heatmap.png` | OG クラスター別運命確率ヒートマップ |

---

### Step 7: 局所 velocity 可視化（ループ検証）

velocity のループが真の生物学的現象か解析上のアーティファクトかを検証。

```bash
python run_og_local_velocity.py
```

**Velocity confidence（局所一貫性スコア）:**

| クラスター | 平均 | 解釈 |
|:-:|:-:|---|
| OG1 | 0.868 | 高い（一方向性） |
| OG2 | 0.865 | 高い |
| OG3 | 0.801 | 高い |
| OG4 | 0.768 | 高い |

全クラスターで 0.75 以上（ループ・ノイズの閾値 0.4 を大きく上回る）。局所的には各細胞の velocity は一貫した方向を向いており、UMAP 投影による見かけ上のループと考えられる。

**主な出力:**
| ファイル | 内容 |
|---|---|
| `figures_og_exact/local_grid_velocity.png` | グリッド矢印（密度比較） |
| `figures_og_exact/local_smooth_compare.png` | 平滑化パラメータ比較 |
| `figures_og_exact/local_phase_portraits.png` | 位相ポートレート |
| `figures_og_exact/pca_stream.png` | PCA 空間での velocity |
| `figures_og_exact/local_velocity_confidence.png` | coherence スコア |

### Step 8: クラスタサイクル検出

OG1–OG4 + PO1–PO2 のクラスタ間遷移から有向サイクルを検出し、サイクルを構成するエッジを赤色で強調表示します。

```bash
python run_og_cycle_detection.py
```

**主な出力:**
| ファイル | 内容 |
|---|---|
| `cycle_detection/cluster_cycle_strengths.png` | クラスタサイクルとエッジ強度の可視化 |

### Step 9: cell-cycle scoring の検証

OG1–OG4 + PO1–PO2 の cell-cycle score を計算し、UMAP、ボックスプロット、クラスター平均値グラフを作成するスクリプトです。

```bash
python run_cell_cycle_scoring.py
```

**注意:**
- 現在の `E17_og_exact_velocity.h5ad` には canonical cell cycle marker genes が含まれておらず、`run_cell_cycle_scoring.py` は有効な遺伝子が見つからない場合に停止します。
- エラーが発生した場合は、`figures_cell_cycle/cell_cycle_error.txt` に詳細を出力します。

---

### Step 10: Loom ファイル統合とバッチ補正（メモリ効率化）

3 つの Loom ファイル（`/home/user/share/SutureRNAVelocity/loom_output/` に格納）を統合し、バッチ補正を行った上で OG1–4, PO1–2 クラスタに絞り込む解析パイプラインを実装。

```bash
# Step 10a: Loomファイル統合 + バッチ補正
python integrate_loom_harmony_og_filtered.py

# Step 10b: Velocityグラフ計算（Dynamicalモデル）  
python compute_velocity_robust.py

# Step 10c: ストリームプロット・サイクル検出・統計分析
python run_final_analysis.py
python generate_streamplot.py
```

**処理フロー:**

1. **Loom 統合**: 3 ファイル（batch_1, 2, 3）を `ad.concat()` で結合（27,849 cells × 55,401 genes）
2. **メタデータマッピング**: `E17_og_exact_velocity.h5ad` の OG クラスター情報を利用して、バーコード照合で 2,017 細胞を抽出
3. **フィルタリング**: OG1–4 + PO1–2 のみを保持（1,743 cells/after QC）
4. **メモリ削減**: Highly variable genes（2,858 遺伝子）のみを保持
5. **バッチ補正**: BBKNN を用いて batch_2, batch_3 のバッチ効果を除去
6. **Velocity 計算**: Dynamical モデル（recover_dynamics）で安定的に推定
7. **ストリームプロット**: streamlines で velocity フローを可視化
8. **サイクル検出**: CellRank VelocityKernel の遷移行列から NetworkX で有向サイクルを検出
9. **統計・可視化**: クラスタ組成、バッチ分布、UMAP+velocity overlay を出力

**主な出力:**

| ファイル | 内容 |
|---|---|
| `E17_og_integrated_harmony_filtered_velocity.h5ad` | バッチ補正済み Loom 統合データ |
| `figures_integrated_analysis/scvelo__stream.png` | Streamlines plot（velocity flow） |
| `figures_integrated_analysis/cycle_detection.png` | クラスタサイクル（20個検出）|
| `figures_integrated_analysis/velocity_overlay.png` | UMAP + quiver plot |
| `figures_integrated_analysis/cluster_composition.png` | クラスタ細胞数分布 |
| `figures_integrated_analysis/batch_distribution.png` | バッチ×クラスタのクロステーブル |
| `figures_integrated_analysis/summary_statistics.csv` | 統計サマリー |

**成果:**

- メモリ使用量を **64GB → < 16GB** に削減（HVG + クラスタフィルタ）
- 3 バッチからの Loom 統合を BBKNN でバッチ補正
- Dynamical モデルで Velocity グラフ計算を成功させた
- クラスタ間の 20 個の有向サイクルを検出（OG2, OG3, OG4 が多く関与）
- バッチ効果の解析が可能に（batch_2: 534 cells, batch_3: 1,209 cells）

---

## ファイル構成

```
SutureRNAVelocity/
├── README.md                          # 本ファイル
├── RNA_velocity_workflow.md           # 詳細ワークフロー
├── aria2c_download.txt                # BAM ダウンロードリスト
├── run_velocyto.sh                    # velocyto 実行スクリプト
├── run_scvelo.py                      # 全細胞 scVelo 解析
├── run_og_exact.py                    # OG1-4 サブセット解析（確定版）
├── run_og_cellrank.py                 # CellRank 脱分化解析
├── run_og_local_velocity.py           # 局所 velocity 可視化
│
├── gencode.vM25.nochr.gtf             # 解析用 GTF（chr プレフィックス除去済み）
├── E17_barcodes.tsv                   # E17 バーコードリスト
│
├── GSE163693_RAW/                     # GEO 発現行列（10x v2 形式）
├── GSE163693_E15_17_composite_metadata.csv.gz  # 論文クラスターラベル
├── downloads/                         # BAM ファイル
├── loom_output/                       # velocyto 出力 loom ファイル
│
├── E17_suture_velocity.h5ad           # 全細胞解析済み AnnData
├── E17_og_exact_velocity.h5ad         # OG1-4 サブセット AnnData
├── velocity_driver_genes.csv          # 全細胞ドライバー遺伝子
├── og_exact_velocity_driver_genes.csv # OG ドライバー遺伝子
│
├── figures/                           # 全細胞解析の図
└── figures_og_exact/                  # OG1-4 サブセット解析の図
```

---

## 主な知見

1. **E17 頭蓋縫合の骨形成軌跡**: OG1（前駆細胞）→ OG2/OG3（プレ骨芽細胞）→ OG4（成熟骨芽細胞）という分化軌跡を RNA velocity で確認。

2. **分岐構造**: OG1 から縫合部（OG2）と骨膜（OG3）の 2 経路への分岐が PAGA で確認された（OG1-OG2 間: 0.772、OG1-OG3 間: 0.752）。

3. **脱分化シグナルの存在**: CellRank GPCCA により OG1 が terminal state の一つとして同定。OG4 細胞の 73% が OG1 方向への velocity を持ち、成熟骨芽細胞から前駆細胞への脱分化経路が示唆される。

4. **velocity の一貫性**: 全 OG クラスターで velocity confidence ≥ 0.77 と高く、局所的な velocity は一方向性を示す。UMAP 上のループ状パターンは投影アーティファクトの可能性が高い。
