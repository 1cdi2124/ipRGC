% =================================================================
% X用(変換実行後)
% =================================================================

[file_name, file_path] = uigetfile('*.set', '常時点灯の.setファイルを選択してください');
EEG = pop_loadset('filename', file_name, 'filepath', file_path);

% パラメータ設定
interval_sec = 0.6; % 周期を他の実験条件（1.2s）に合わせる
srate = EEG.srate;
interval_pnts = interval_sec * srate;
total_pnts = EEG.pnts;

% イベントの生成（データの端に余裕を持たせる）
EEG.event = []; % 既存のイベント（もしあれば）をクリア
count = 1;
for p = interval_pnts : interval_pnts : (total_pnts - interval_pnts)
    EEG.event(count).type = 'Saccade'; % 解析スクリプトが探す名前 'Saccade' に合わせる
    EEG.event(count).latency = p;
    EEG.event(count).duration = 0;
    count = count + 1;
end

EEG = eeg_checkset(EEG);
% 保存（元のファイル名に _dummy と付けて保存）
[~, base, ~] = fileparts(file_name);
pop_saveset(EEG, 'filename', [base '_events.set'], 'filepath', file_path);
fprintf('イベント付与完了: %d 個のイベントを設置した。\n', count-1);