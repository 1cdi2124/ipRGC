# ERSP品質スクリーニング・試作版

`run_ersp_quality_screening.m` は既存の共通ICA/newtimef出力を読み、品質不良の**候補**を抽出し、既存06帯域出力から条件×帯域×チャンネルの三要因ANOVA用CSVを作る独立プログラムです。教授による画像評価の補助として、品質閾値はまだ未検証です。`CWT_20261006` は04～06を作る役割に分離しました。

## 実行

MATLABのCurrent Folderを、このプログラムを置いたフォルダにしてください。EEGLABは起動不要です。

```matlab
run_ersp_quality_screening
```

フォルダ選択画面で次のいずれかを選びます。

- 被験者一覧（例：`01_CWT_joint_new_NoGreen_Brain`）：直下の被験者全員。
- 被験者フォルダ：その被験者のみ。
- `03_手動サッケード` / `04_ERSP_共通ICA_newtimef_全帯域` / `06_バンド帯_共通ICA_newtimef_全帯域` / その中の電極フォルダ：その被験者の指定条件・指定電極。

異なる実験や旧解析版を混ぜないため、深い再帰探索はしません。複数の解析版を含む上位の実験フォルダではなく、**1つの解析版の被験者一覧**を選んでください。`hasegawa_0807_F` と `hasegawa_0820` はフォルダ名をそのまま別IDとして扱います。

```matlab
% 保存先を別の親フォルダへ指定
[cells, subjects, out] = run_ersp_quality_screening("C:\...\被験者一覧", "D:\QC")

% 色範囲・試行数警告の仮閾値を変える例
options = struct('display_limit_db', 1.5, 'min_trials_review', 30);
run_ersp_quality_screening("C:\...\被験者一覧", "", options)

% 3条件に限定する場合（初版の既定は4条件）
options = struct('conditions', ["0 Hz", "80 Hz", "160 Hz"]);
run_ersp_quality_screening("C:\...\被験者一覧", "", options)
```

第1引数：入力。省略・空文字で選択画面。第2引数：出力の親フォルダ。省略・空文字で対象範囲の `08_ERSP_QC`。第3引数：設定を上書きするscalar struct。設定の全既定値はコード内の `default_config()` にあります。元の解析プログラムのconfigは変更しません。

## 実行モードと07_ANOVA（2026-10-07追加）

`options.mode` は次の4つです。

- `"all"`（既定）：04の品質評価・HTML画像一覧と、06から07の三要因表作成。
- `"qc_only"`：品質評価・HTML画像一覧だけ。07を作らない。
- `"anova_only"`：既存06から07だけを作る。ICA、ERSP計算、QC画像生成を繰り返さない。
- `"wide_only"`：`ANOVA_all_cells_QC.csv` のあるフォルダを選択し、既存集計値から三要因表だけを整形。旧CWTのwide_onlyはこちらへ移動。

```matlab
% 全員のQCと三要因表（被験者一覧フォルダを選択）
run_ersp_quality_screening

% ANOVA用CSVだけを作る（画像生成の待ち時間を省く）
run_ersp_quality_screening("", "", struct('mode', "anova_only"))

% 以前と同じQCだけ
run_ersp_quality_screening("", "", struct('mode', "qc_only"))

% 既存cells表だけを整形
run_ersp_quality_screening("C:\...\cells表のあるフォルダ", "", struct('mode', "wide_only"))
```

07の既定先は、選択範囲の `07_ANOVA_バンド帯_チャンネル/run_日時/` です。第2引数で親フォルダを指定した場合は、その親の下に07を作ります。`options.anova_output_root` で07の親を別途指定することもできます。どのモードでも既存のANOVA表や手動判定は上書きしません。`anova_only` / `wide_only` の第3戻り値 `out` は07の保存先で、最初の2戻り値は空です。

時間窓は `config.anova_time_window_ms = [0 300]`。4条件（0/80/160 Hz/NoGreen）×8帯域×8電極が既定です。`conditions` / `channels` / `anova_band_names` で因子を指定できます。帯域名や列順は従来の三要因表を維持しています。**条件別の二要因表は作りません。**

06 MATの `band_mean` / `band_n` / `times` を優先し、同じstemのCSVは二重計上しません。MATがないファイルのみ06 CSVの `MeanERSPdB` / `TrialCount` を使います。既存の帯域dB値を時間窓内で平均し、dB変換や周波数平均を再実行しません。MATが壊れている場合、同じstemのCSVへ黙って切り替えずinvalidとして記録します。

06がない条件は、04の保存済み `ersp_mean` と元の `analysis_config.band_ranges` を使い、CWTの `align_band_means_with_newtimef` と同じ周波数方向のdB平均を作ります。ICA・newtimef・dB変換・補間は再実行しません。指定時間窓の該当周波数で `ersp_n == trial_count` がすべて成立する場合のみ、有効試行数を正確に復元します。一部試行の欠測対応が不明な場合は推測で埋めずinvalidとし、06が必要である旨を記録します。読取不能な06を04へ黙って置き換えません。どちらを使用したかはMATのSourceFilesに記録します。

07のCSVは次の5表です。

- `ANOVA_all_cells_QC.csv`：全員の条件×帯域×chの集計値、時間点数、最小試行数、時間窓。
- `ANOVA_threeway_wide.csv`：SPSS用。1行1被験者、先頭ParticipantID＋因子順に並ぶ測定列（既定は256列）。独立したCondition列は付けません。
- `ANOVA_threeway_factor_design.csv`：測定列と条件・帯域・chの対応、SPSSへ登録する順番。
- `ANOVA_threeway_inclusion_QC.csv`：完全データかどうか、欠測・重複・invalid、保存データから取得したICA設定と参考用のQC候補。
- `ANOVA_threeway_long_complete.csv`：wideへ採用された完全データを縦形式で保持。

`ANOVA_export_results.mat` に読み込んだ04/06の参照や設定を記録します。SourceCSVやOriginalCondition/SessionIDの重複列はANOVA表へ追加しません。CSVのみの場合、ICA設定はUnknown/NaNとし、現在のCWT設定で解析済みだとは扱いません。手動ICAの場合、Brain閾値は初期候補の基準であり、自動採用閾値としては記録しません。

**QCの要確認・除外候補で被験者を自動除外しません。** wideへの採用は従来どおり、指定因子の全セルが一意で値・試行数・時間窓が有効な被験者です。QC候補はinclusion表へ参考として記録します。確定した品質基準での除外は別途判断してください。単被験者のwideは作れますが、群の反復測定ANOVAには複数被験者が必要です。

## 出力

毎回 `08_ERSP_QC/run_年月日_時分秒_ミリ秒/` を新規作成します。前回結果や手動判定を上書きしません。

- `QC_cells.csv`：被験者×条件×電極の指標、候補、理由、入力MATの参照。
- `QC_participants.csv`：被験者別の警告・欠測セル数と候補。`FinalDecision` / `Reviewer` / `FinalReason` は手動判定用の空欄。
- `QC_review.html`：採用候補も含む全セルの画像一覧。既定で「全て」を表示し、「全て／両方／要確認／除外候補」の4つから切り替えられます。「両方」は要確認＋除外候補です。被験者・条件・電極で検索できます。画像をクリックして拡大できます。
- `images/`：採用候補も含む全セルの再描画画像。全セルに共通の色範囲・軸で、補間せず描きます。警告セルだけを描く旧版より、画像作成に時間がかかります。
- `boards/`：警告被験者について、**警告なしも含む**全条件×全電極の一覧画像。ファイル名は被験者フォルダ名そのままの `fujiyama_0827.png`、`hasegawa_0807_F.png`、`hasegawa_0820.png` などです。日付・接尾辞は削らず、別人を同じ画像へ統合しません。HTMLとCSV/MAT内の一覧画像参照もこの名前を使用します。Missing/duplicateも明示します。
- `QC_results.mat`：指標表、入力フォルダ、実行時設定、実行開始時のプログラムソースと07出力情報を記録。処理の追跡用です。

CSVがExcelの通常ダブルクリックで文字化けする場合は、UTF-8を指定して取り込んでください。HTMLはUTF-8を明示しています。HTMLを閲覧してもネットワークには接続しません。

旧版で作成済みのHTMLは自動では変わりません。更新したプログラムを再実行して、新しい `run_日時/QC_review.html` を開いてください。ICAの再計算は不要です。品質指標・採用候補の判定条件は、今回の全件表示への変更では変えていません。

## 指標と判定の意味

判定に用いる形状は、元MATの `analysis_config.baseline_window_ms` に対応する**イベント前だけ**です。元MATに設定がなければ `[-250 -100] ms` を仮使用し、要確認になります。周波数は既定で30–200 Hz。低周波の時間的なにじみを、縦筋の異常判定と混同しないためです。0.5–30 Hzも画像には表示しますが、この初版では形状判定が不十分です。

| 指標 | 定義 | 仮の要確認閾値 | 仮の強い異常閾値 |
| --- | --- | --- | --- |
| BaselineVerticalDb | 中心列と左右10–30 msの列の中央値平均との差の絶対値を、周波数方向に中央値化した値の最大 | 0.75 dB | 2.0 dB |
| BaselineHorizontalDb | 中心行と上下3–10 Hzの行の中央値平均との差の絶対値を、時間方向に中央値化した値の最大 | 0.75 dB | 2.0 dB |
| BaselineRoughnessDbPer10Ms | 隣接時間点の絶対dB差を時間差で割り、全点の中央値を10 ms当たりに換算 | 1.0 dB/10 ms | 3.0 dB/10 ms |
| ChannelRoughnessRatio / ChannelRobustZ | 同じ被験者・同じ条件の電極間で粗さを比較。比とlog10指標のrobust Z | 比3以上かつZ3.5以上 | 単独では除外候補にしない |

これらは**試作の形状指標であり、既存の確立したEEGノイズ尺度ではありません**。指標の「粗さ」は平均ERSPの不規則性で、電極抵抗やEEGノイズ量の直接推定ではありません。教授が選んだ良い・悪い・判断が難しい例に照合して、閾値を検証してください。指標は相関し得るため、2項目が高いことを独立した2つの証拠とみなせません。

- `AcceptCandidate`＝採用候補：本指標で警告がない。品質保証や確定採用ではありません。
- `Review`＝要確認：仮閾値の超過、20試行未満、baseline情報不足、電極間の粗さの突出など。
- `ExcludeCandidate`＝除外候補：非有限値、不正・読取不能MAT、2試行未満、重複セル、または複数の強いbaseline形状異常。**自動除外はしません。**

被験者は、期待する条件×電極のセルのうち除外候補が25%以上なら除外候補、それ以外の警告・欠測があれば要確認になります。割合25%、試行数20などは仮設定で、研究に合わせて固定する必要があります。1つの警告セルだけで被験者全体を除外候補にする設計ではありません。欠測とノイズは区別して記録します。

`DisplayVerticalDb` / `DisplayHorizontalDb` は共通表示区間全体の形状の参考値です。イベント後の本来の活動で高くなる可能性があるため、候補判定には使いません。`DisplayClippedFraction` は色範囲を超える割合で、**品質判定には使いません**。色範囲を変更しても指標・候補は変わりません。

## 制限と安全性

- 入力は `*_共通ICA_newtimef_ERSP.mat`。`ersp_mean`（周波数×時間）、`times`、`freqs` が必須。PNGしか残っていない場合は非対応です。
- 試行平均ERSPだけでは、少数の不良試行、飽和、持続する電源ノイズ、電極抵抗、イベント誤検出を確定できません。疑わしいものは元EEG・PSD・試行数等で確認します。
- 元ERSPがbaseline補正された相対量なので、持続するノイズが見えないことがあります。
- 時間周波数変換の時間的なにじみによって、イベント後の活動がイベント前の値へ影響する可能性もあります。baselineだけで判定しても、生理反応との完全な分離は保証されません。
- 波形、SET、ICA採用条件、解析値、既存ANOVA用wideを変更しません。新規07に表を作りますが、データの削除・試行除外・補間・ICAの再計算はありません。
- 本来の反応や個人差を「不良」と誤判定し得るので、色の強さ・条件差・p値で閾値を調整しないでください。
- 今回の統計結果を見た後に導入したQCであることを記録し、17名の元の結果と品質基準適用後の結果を保存してください。

## 検証

2026-10-07に `01_CWT_joint_new_NoGreen_Brain` の実データを読取りのみで検証しました。17名×256測定列を作成でき、現在の06 MATの指定時間窓平均と一致しました（CSVとの差も丸め誤差の範囲）。06のない `takahira_0807` は04に保存された値から同じ帯域平均を復元できました。既存の集計CSVを `wide_only` で整形した表は、旧wideと完全一致しました。ただし、旧集計値と現在の04/06から作る値は一致しません。現在データの再集計と旧表の整形は別の操作なので、旧SPSS結果に新wideをそのまま対応づけないでください。旧表・元データは変更していません。

ICA前後の保存警告は、EEGLAB `pop_saveset` にstring型の保存先を渡すと内部のファイル名が配列になり、短絡論理演算でエラーになるためでした。CWTの保存先をchar型へ変換して修正しています。ICA前後の出力を作り直すには更新したCWTの再実行が必要です。ICA採用基準・seed・ERSP計算は今回変更していません。

07移動・ICA保存修正のテストは `C:\Users\tatsuya\Documents\ipRGC\cwt_qc_anova_20261007` にあります。`test_qc_threeway_integration` はMAT/CSV混在、欠測、重複、単被験者、既存表整形、QC除外候補を自動除外しないことを確認します。`test_cwt_qc_output` は小さな合成EEGでICA前後のSET/FDT・波形・境界をまたがないPSDを検証します。研究データのICAを再実行するテストではありません。

合成テスト `test_ersp_quality_screening.m` は開発用作業フォルダ `C:\Users\tatsuya\Documents\ipRGC\ersp_qc_prototype_20261006` に置いています。テスト専用フォルダに合成MATと出力を作り、既存の研究データを変更しません。共有側の解析実行コードには本体とこの使い方だけを配置しています。

```matlab
addpath('C:\Users\tatsuya\Documents\ipRGC\ersp_qc_prototype_20261006')
test_ersp_quality_screening
```

非有限値、少数試行、縦横筋、欠測・重複、別被験者ID、単一被験者選択、採用候補を含む全セルのHTML/PNG出力、4つの切替選択肢、強いイベント後反応だけで除外しないことを確認します。これはプログラム動作の検証であり、品質分類の科学的妥当性の検証ではありません。

参考：EEGLABの視覚的確認との併用という考え方と、その旧版手法の説明：https://eeglab.org/tutorials/misc/Rejecting_Artifacts_Legacy_Menus.html （この試作コードは旧版の拒否関数を呼びません）。
