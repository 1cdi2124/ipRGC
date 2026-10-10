# ipRGC

Phantom Array Effect（PAE）研究の解析コードを管理するリポジトリです。
「04　解析実行コード」フォルダーをGitHubの既存履歴に接続しています。

## コードの説明

[解析コード用途ガイド](00_解析コード用途ガイド/README.md)と[詳細説明](00_解析コード用途ガイド/説明.md)を参照してください。
既存の `dummy_event.m` は履歴との互換性のため残しています。共有フォルダー由来の `dummy_event_.m` とは別ファイルです。
実行環境やデータのパスは各コード内の設定を確認してください。

## 変更を保存する手順

このフォルダーでPowerShellを開き、以下を実行します。

```powershell
# 作業前。他のPCなどが更新した履歴を取得（未保存の編集がある場合は先にcommit）
git pull --ff-only

# MATLAB等で編集後、差分と対象を確認
git status
git diff

# 変更したファイルを指定して記録・送信（ファイル名は例）
git add -- run_joint_ica_newtimef_band_analysis.m
git diff --cached
git commit -m "解析設定を更新"
git push

# 変更履歴を確認
git log --oneline -20
```

ファイルの保存やOneDriveの同期だけではGitHubに履歴は登録されません。変更の区切りごとにcommitとpushを行ってください。
同じ共有フォルダーを複数PCで同時にGit操作しないでください。別PCではOneDrive外へ個別にcloneして利用する構成を推奨します。

## 管理対象

このフォルダー内のコード・説明文書・データ・ZIPを管理します。2026年9月26日にデータとZIPを含む公開の承認を受けて登録しました。
`M/tokai-internship-main.zip` はデータや資料を含む元のアーカイブのまま保存します。ZIP内の個別ファイルの差分はGitHub上では確認できません。
OSの管理ファイルやPythonのキャッシュ、仮想環境などは `.gitignore` で除外します。
このリポジトリは公開されています。新たに追加するファイルは、commit前に `git status` と内容を確認してください。
