## 2026-03-05 セッション履歴（Cursor）

- **OG サブセットの拡張**
  - もともと OG1–4 のみを対象としていたサブセット解析 (`run_og_exact.py`, `run_og_cellrank.py`, `run_og_local_velocity.py`) を、**PO1 / PO2 を含む 6 クラスター（OG1–4 + PO1–2）** に拡張した。
  - メタデータ `GSE163693_E15_17_composite_metadata.csv.gz` から E17hiseq の OG/PO バーコードを再取得し、`E17_suture_velocity.h5ad` から 2,017 細胞をサブセット化。

- **CellRank 解析のループ許容化とラベル統一**
  - `run_og_cellrank.py` で VelocityKernel + ConnectivityKernel を用いて CellRank GPCCA を実行し、**terminal/initial state による一方向の吸収モデルではなく、ループを含むトポロジー解釈**に変更。
  - macrostates（OG3_1, OG3_2, OG4_1 など）は内部では保持しつつ、プロットでは常に **元のクラスター名 `og_cluster`（OG1–4, PO1–2）で色付け**するように変更（`macro_major` に多数派ラベルをマッピング）。
  - PAGA 接続行列と、クラスター間の正味フラックスに基づく **粗視化された分化ベクトル場 (`cr_diff_vector_field.png`)** を作成。

- **局所 velocity と信頼度の 6 クラスター対応**
  - `run_og_local_velocity.py` を OG1–4 + PO1–2 対応に修正し、各クラスターを個別にハイライトした局所ストリーム、PCA 空間での velocity、velocity confidence のボックスプロットなどを再計算。

- **PowerPoint 自動生成の更新**
  - 解析サマリー用スライド (`make_pptx.py`) を、**PO1 / PO2 を含む 6 クラスター前提**に更新。
    - クラスター定義スライドを「OG/PO Cluster Definition」に変更し、PO1 / PO2 の行を追加。
    - OG/PO サブセット UMAP・velocity・latent time・velocity confidence の説明文・値を 6 クラスター版に更新。
  - `figures_og_exact/` 配下の全図（OG/PO サブセット解析）を 1 枚ずつ並べる専用の PowerPoint を作るスクリプト `make_pptx_og_exact.py` を新規追加。

- **Mmp13 発現の再取得と可視化**
  - `E17_og_exact_velocity.h5ad` には Mmp13 が含まれていなかったため、**元の 10x 行列 `GSE163693_RAW/`（`GSM4983999_E17_*`）から Mmp13 カウントを再取得**。
  - OG/PO サブセットとのバーコードを照合し、各細胞ごとの `Mmp13_expr` と、その 5–95 パーセンタイルでクリップして 0–1 にスケーリングした `Mmp13_scaled` を `obs` に追加。
  - UMAP + velocity の配置は `run_og_exact.py` と同じものを用い、**Mmp13 発現（スケーリング後）で色付けしたストリーム図** `scvelo__og_exact_Mmp13_scaled.png` を `figures_og_exact/` に出力。

- **figures_og_exact の再生成**
  - 途中で図が増えすぎたため、`figures_og_exact` 配下をいったん削除し、修正済みの `run_og_exact.py` を再実行して、OG/PO サブセットの図一式（UMAP、velocity、latent time、proportions など）と Mmp13 図を再生成した。

---

## 2026-03-25 セッション履歴（Copilot）

- **スライド 8-9 の PAGA・CellRank 解析についての理論的解説**
  - PAGA グラフの構造：無向グラフでクラスター間の接続を表す。エッジの太さ = クラスター間の遷移強度。
  - CellRank の役割：RNA velocity と PAGA 接続情報を統合してクラスター間の運命遷移を推定。
  - Initial / Terminal state の決定方法：CellRank GPCCA が定常分布を計算し、自動的に識別。このコードでは Terminal: OG4 + OG1、Initial: OG3。
  - 遷移確率行列の構築：VelocityKernel（80%）+ ConnectivityKernel（20%）で加重統合。
  - PAGA 接続の算出方法：各細胞の k-nearest neighbor グラフから、クラスター間の統計的に有意な相互作用（流出）を検出。

- **Initial/Terminal state が生物学的知見と不一致する原因の検討**
  - OG1 が初期かつ終末状態として識別される問題を指摘。
  - 考えられる原因：
    1. Kernel 加重比率（VelocityKernel 80% ← 高すぎる可能性）
    2. scVelo 動力学モデルの収束不足（特に小規模クラスタ）
    3. クラスタサイズの不均衡による bias
    4. **Batch effect の存在**（3 バッチ SRR13288596–98 からのデータ）
    5. GPCCA macrostate 数の設定不適切
    6. Temperature パラメータの最適化不足

- **Batch effect の説明と補正**
  - Batch effect：異なる実験条件下（異なる時間、プレート、試薬ロット）で生成されたデータ間の技術的系統誤差。
  - このコードの場合：3 つの SRR run（3 バッチ）からのデータが batch effect を含む可能性。
  - Batch effect の症状：RNA 発現量、RNA velocity 方向が逆向きになる等。
  - PAGA 接続性や velocity の推定精度に悪影響。

- **Batch correction 対応への修正**
  - `run_og_cellrank.py` を以下のように修正：
    1. Harmony または BBKNN による batch correction の追加
    2. Batch 情報の自動検出（`batch` 列および `SRR` / `run_accession` 列）
    3. Batch correction 前後の UMAP 比較図の出力
    4. 新しい出力ディレクトリ `figures_og_exact_batch_corrected/` に結果を保存
  - 依存パッケージを `pyproject.toml` に追加：`bbknn`, `harmonyphpy`

- **環境構築と実行準備**
  - Python 仮想環境 `.venv` を新規作成
  - 依存パッケージをインストール：numpy, scipy, scanpy, scvelo, cellrank, bbknn, harmonypy
  - `run_og_cellrank.py` の batch correction 対応版の実行準備完了
  - 次回：実際に batch correction 解析を実行し、OG1 が terminal state として識別されるかどうかを検証予定

- **2026-04-04 セッション追加**
  - `run_og_cycle_detection.py` を追加し、OG1–OG4 + PO1–PO2 のクラスタ間遷移で有向サイクルを検出・可視化する機能を実装。
  - サイクルは NetworkX の `simple_cycles()` で検出し、サイクル内エッジを赤色で強調表示する可視化を `cycle_detection/cluster_cycle_strengths.png` に保存。
  - `run_cell_cycle_scoring.py` を追加し、OG1–OG4 + PO1–PO2 の cell-cycle score を算出して UMAP、ボックスプロット、クラスター平均値グラフを出力する機能を追加。
  - ただし、`E17_og_exact_velocity.h5ad` には canonical cell cycle marker genes が含まれておらず、現状では cell-cycle scoring が計算できないことを確認。
  - スクリプトは有効なマーカー遺伝子が存在しない場合に明示的に停止し、`figures_cell_cycle/cell_cycle_error.txt` にエラーメッセージを書き出すように修正。

---

## 2026-04-10 セッション履歴（Copilot）

### 背景
- サイクル検出スクリプト (`run_og_cycle_detection.py`) でOG1とOG3の循環遷移が検出されたが、ストリームプロット (`scvelo__og_exact_stream.png`) には当該経路が明確に見えないという不一致が報告された。
- 原因の根本調査と、バッチ補正付きの Loom 統合パイプラインの構築を実施。

### 1. サイクル検出の不一致原因分析
  - **検出結果**: CellRank VelocityKernel + ConnectivityKernel（加重: 0.8 + 0.2）で OG1 ↔ OG3 の双方向遷移を検出
  - **ストリームプロット**: scVelo の velocity embedding stream では OG3→OG1 経路が視覚的に目立たない
  - **結論**: 
    1. ConnectivityKernel（k-NN近傍グラフ）の影響が大きい可能性
    2. Velocity ベクトルの強度が弱い（流線として目立たない）
    3. UMAP 投影での非線形変換によるマスキング
    4. **バッチ効果の存在**（3 バッチからの混在データ）

### 2. Loom ファイル統合スクリプトの作成
  - **構築スクリプト**: `integrate_loom_harmony_og_filtered.py`
  - **特徴**:
    - `/home/user/share/SutureRNAVelocity/loom_output/` の 3 Loom ファイルを `ad.concat()`で統合（27,849 cells × 55,401 genes）
    - メモリ削減のため **HVG（高変動遺伝子）2,858 個のみを保持**
    - バッチキー `batch` を自動検出（batch_1, batch_2, batch_3）
    - 既存の `E17_og_exact_velocity.h5ad` のメタデータ（og_cluster）を利用してバーコード照合
    - **OG1–4, PO1–2 クラスタ（1,743 細胞）のみをフィルタリング**
    - BBKNN によるバッチ補正（Harmony インストール試行後、API 非互換で BBKNN に変更）
    - 出力: `E17_og_integrated_harmony_filtered.h5ad`

### 3. Velocity 計算の安定化
  - **問題**: stochastic モデルで "setting an array element with a sequence" エラー
  - **原因**: スパースデータ（90% 以上がゼロ）で最小二乗法が数値的に不安定
  - **解決策**: Dynamical モデル（`recover_dynamics` + `velocity(mode="dynamical")`）を採用
  - **構築スクリプト**: `compute_velocity_robust.py`
  - **結果**: Velocity 計算成功、velocity_graph と velocity レイヤが正常生成

### 4. メモリ削減のポイント
  - **元の課題**: 27,849 cells × 55,401 genes で 64GB メモリ使用、Harmony 実行時に強制終了
  - **削減戦略**:
    1. **クラスタフィルタ**: OG1–4, PO1–2 のみ（27,849 → 2,017 細胞）
    2. **HVG 選択**: min_mean=0.0125, max_mean=3, min_disp=0.5（55,401 → 2,858 遺伝子）
    3. **結果**: < 16GB メモリで処理完了

### 5. ストリームプロット + サイクル検出 + 統計分析
  - **ストリームプロット**: 
    - scVelo `velocity_embedding_stream()` は失敗（'dict' object has no attribute 'dtype'）
    - 代替実装: scipy `griddata` で velocity field を補間し、matplotlib `streamplot` で可視化
    - 出力: `figures_integrated_analysis/scvelo__stream.png`
  - **サイクル検出**: CellRank VelocityKernel で遷移行列を構築、NetworkX で有向サイクルを検出
    - **検出結果**: **20 個のサイクルを検出**
    - 主要サイクル: OG2 ↔ OG1 (0.212/0.174), OG2 ↔ OG3 (0.192/0.160), OG4 を含む複雑な経路
    - 出力: `figures_integrated_analysis/cycle_detection.png`（赤色でサイクル内エッジを強調）
  - **統計分析**:
    - クラスタ組成: OG3 (527), OG1 (399), OG2 (288), OG4 (288), PO1 (202), PO2 (39)
    - バッチ分布: batch_2 (534), batch_3 (1,209)
    - 出力: cluster_composition.png, batch_distribution.png, summary_statistics.csv

### 6. 最終出力ファイル
  - **統合データ**: `E17_og_integrated_harmony_filtered_velocity.h5ad` (311 MB)
  - **可視化**:
    - `scvelo__stream.png`: ストリームプロット
    - `cycle_detection.png`: クラスタサイクル（20 個）
    - `velocity_overlay.png`: UMAP + quiver plot
    - `cluster_composition.png`: 棒グラフ
    - `batch_distribution.png`: 積み上げバー
  - **メタデータ**: summary_statistics.csv, batch_cluster_distribution.csv, cluster_composition.csv

### 7. 主な学習ポイント
  - **バッチ補正の段階**: Loom 統合後、前処理前に適用することが重要（メモリ効率と精度のバランス）
  - **Velocity モデル選択**: scVelo のデフォルト stochastic はスパースデータで不安定。Dynamical モデル推奨
  - **scVelo API の制限**: velocity_embedding_stream に不具合の可能性。手動実装で回避可能
  - **HVG 保持**: 生物学的に有意な遺伝子のみを保つことで、計算量削減と解釈性向上が両立
  - **メモリ効率化**: クラスタ + 遺伝子のダブル制約で 64GB → 16GB 以下に削減

### 8. 今後の課題
  - バッチ補正後のサイクルパターンの変化を追跡（artifact vs. 生物学的現象の判定）
  - CellRank GPCCA による macrostate 検出の再実行（バッチ補正後）
  - Harmony インストール（現在 harmonyphpy が pip で取得できるが、Scanpy 統合 API が未対応）
  - より詳細な cell-cycle marker gene の拡張（外部データベース統合）

