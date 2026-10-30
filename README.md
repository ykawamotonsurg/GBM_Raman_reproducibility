# GBM Raman: minimal saved-data plotting package

**最終図の描画に必要な保存データとコードの最小公開セットです。S2Bの採用結果と描画コードを含みます。**

最終図の描画に必要な保存済みデータとコードをまとめています。データCSVは30個です。完成PNG、画像配置用の素材、内部監査資料、原稿、患者実日付・対応表は同梱しません。補足表は原則として別提出の `Data S1.xlsx` を参照し、描画または保存された注記に必要な表だけCSVとして含めます。

## 主図から探す

| 図 | 内容 | 主なデータ |
|---|---|---|
| Fig.1 | 組織スコア・ROC・スペクトル・ピーク | `main/data/tissue_bags.csv`, `main/data/tissue_spectra.csv`, `shared/data/tissue_*.csv` |
| Fig.2 | 単一細胞スペクトル・保存UMAP・high割合 | `main/data/cell_mean_spectra.csv`, `shared/data/cells.csv`, `main/data/Singlecell_CaseRegion.csv` |
| Fig.3B–D | 同一細胞対応・培養実験 | `main/data/Table_S6.csv`, `Culture_AllConditions.csv`, `Fucoidan_Pairs.csv`, `Table_S7A.csv`, `Table_S7B.csv`（いずれも `main/data/`） |
| Fig.4B | 空間プログラムの保存スコア | `shared/data/spatial_program_scores.csv` |

**各パネル・補足表の正確な対応は [figure_data_map.csv](metadata/figure_data_map.csv) を参照してください。** Fig.3A、Fig.4A、Fig.S1、Fig.S4Aの画像パネルは論文側の図を参照します。このZIPだけで全図の完成レイアウトが復元できるという意味ではありません。

## フォルダ

- `main/data/`：主図用のデータ9個。`main/code/`：主図の描画コード。
- `supplement/data/`：補足図用のデータ15個。`supplement/code/`：補足図の描画コード。
- `shared/data/`：両者が共用するデータ6個。患者対応を保った同じファイルを参照します。
- `metadata/`：図との対応、CSVの列説明、入力ハッシュ、別提出資料の参照先。

[入力ファイル一覧](metadata/input_files.csv) / [列説明](metadata/data_dictionary.csv) / [Data S1との対応](metadata/submission_references.csv)

## 実行

Python 3.12で、ZIPを展開したフォルダから実行します。パッケージの版は `requirements.txt` と `environment.json` に記録しています。

```sh
python -m pip install -r requirements.txt
python reproduce.py --scope all --out replay
```

主図だけなら `--scope main`、補足図だけなら `--scope supplement` を指定します。出力先は新規または空のフォルダとし、再実行時は別名を使います。主図12パネル・補足図15パネルのPNGと、Fig.3Cの別サイズの補助PNGを生成します。`replay/` は生成物なのでGitHubに登録する必要はありません。図の書体・余白等は論文の完成配置と異なる場合があります。

## 再現の範囲

このパッケージは **保存数値からの描画** を行います。生データからの前処理、学習、UMAP計算、閾値選択、統計検定を再実行するものではありません。画像をコピーして原データ再現とみなす処理もありません。

1,499細胞・high 225細胞、主解析57 bags、固定された患者対応・座標・分類・スコア・P値・区間・ICC(2,1)を保持します。箱ひげ、表示用logit、既存の記述的直線等は従来の描画処理を使用します。一部の注記は保存コードに固定された文字列で、根拠の保存結果CSVを併記しています。新しい推定・検定は行いません。

Fig.3Bは `Table_S6.csv` の確定済みカウントを四角で描きます。以前の保存用セットにある `samecell_labels.csv` はこの最小版には含めていません。元の測定からhigh/lowやMSR1判定を作り直せるとは主張しません。Fig.S6Bの上流にある症例間集約の詳細も、この描画用パッケージでは独立に確認していません。

## 別提出資料と補足表の提供先

補足表S1–S9は `Data S1.xlsx` の該当シートを参照します。提出予定の提供先を記したもので、現在その資料が一般公開されているという意味ではありません。公開時には採用された提出資料・リンクと対応表を一致させてください。Data S1を別途提供しない場合は、必要な表を本パッケージへ戻す必要があります。

## S2B：時間分割評価の保存結果

- `supplement/data/s2b_temporal_late_bags.csv`：後期22 bags／16患者の公開用ID、領域、スペクトル数、確定bagスコア、保存分類。
- `supplement/data/s2b_temporal_threshold.csv`：top-3 mean用の固定閾値 **0.98162150382995605** と比較規則。通常表示の0.981621503829956と同じ保存閾値です。
- `supplement/code/plot_s2b.py`：保存分類を領域別に数えてヒートマップを描き、固定閾値との整合を確認します。

```sh
python reproduce.py --scope s2b --out replay_s2b
```

S2BのPNGと領域別カウントCSVを生成します。tumor予測はCE 10/10、IM 2/3、ED 1/3、CWM 0/6です。元コードのNAWMを現行原稿のCWM表記へ揃えたもので、領域判定を変更していません。全16患者のIDを既存の公開用対応から引き継いでいます。

入力は追跡可能な再実行の採用結果です。主解析57 bagsのOOFスコアや旧閾値を混ぜていません。保存記録上、前処理とモデルは前期CE/CWMのみで学習し、閾値は前期IM/EDのOOFのみで決定しています。この公開用コードは学習・推論・閾値決定を実行しません。

旧実行とは22 bagsのtop-3分類が一致しましたが、連続予測値と閾値は一致していません。当時の実行履歴や完全な数値再現を証明するものではありません。ここで再現できるのは保存結果からの描画です。学習済みモデル、前処理、元入力や患者日付の対応表はこの最小セットには含みません。S2全体の完成画像の配置再現も含みません。

今回版は `main/`、`supplement/`、`shared/` と本READMEの対応表を入口にしてください。旧図番号の `fig1_public_data/`、`fig2_public_data/`、`fig3_public_data/` は[更新前のGitHub履歴](https://github.com/ykawamotonsurg/GBM_Raman_reproducibility/tree/1749fdf46532dcf220c8f30cbe5109f10fb738db)から参照できます。旧版と今回版の図番号・スコアを混同しないでください。

今回版を引用する際は、使用したGitHubコミットの固定URLを記録してください。`main` は更新される可能性があるため、特定の保存結果を引用する際はコミットIDで版を区別します。

## 旧DOIと今回版の提供先

既存DOI **[10.5281/zenodo.19363426](https://doi.org/10.5281/zenodo.19363426)** は、旧版 `v1.0.0` の保存先として維持します。このDOIに今回の更新ファイルが保存されているという意味ではありません。今回版はGitHubの確定したコミットを別途参照し、原稿にも旧版DOIと今回版のGitHubリンクを区別して記載します。旧版のアーカイブ・タグ・ReleaseとGit履歴は保持します。この公開方針では新しいZenodo版・DOIを作成しません。

## English scope

This package redraws 27 panels from 30 fixed, pseudonymized CSV inputs. It does not reconstruct the original measurement-to-analysis pipeline. No model training, threshold optimization, UMAP refitting or inferential testing is performed. Completed figure images and image-layout reconstruction are excluded. Supplementary tables are referenced to the separate author-held Data S1 submission; publication links remain to be finalized. Fig. S2B is included as a saved-result redraw of 22 late-period bags from 16 patients, using the adopted top-3 scores, fixed threshold and stored classifications. It does not reproduce upstream training or historical continuous predictions. Use an exact GitHub commit URL to identify the revision used.

The existing DOI, https://doi.org/10.5281/zenodo.19363426, identifies the original v1.0.0 archive. It does not archive the revised files in this package. The revised package is provided through GitHub and must be identified separately by the exact commit used. The existing DOI is retained, and no new Zenodo version is planned under this publication approach.
