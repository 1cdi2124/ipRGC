% =================================================================
% 200ms1sのcsvの解析用
% 緑トリガーの開始地点にイベントを設置するプログラム
% =================================================================

% 1. フォルダの選択
input_dir = uigetdir('', '生データ(CSV)が含まれるフォルダを選択してください');
if input_dir == 0, error('フォルダが選択されなかったため, 処理を中断した.'); end

% 保存先フォルダの作成
[parent_dir, ~, ~] = fileparts(input_dir);
output_dir = fullfile(parent_dir, '02_緑開始点検出');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% グラフ保存用フォルダの作成
plot_save_dir = fullfile(output_dir, 'StimulusChecks');
if ~exist(plot_save_dir, 'dir'), mkdir(plot_save_dir); end

% CSVファイルの取得
files = dir(fullfile(input_dir, '*.csv'));
if isempty(files), error('CSVファイルが存在しない.'); end
fprintf('合計 %d 個のCSVファイルを処理する.\n', length(files));

% EEGLABの初期化
if ~exist('eeglab', 'file')
    error('EEGLABのパスが通っていない.');
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% パラメータ設定
chan_labels = {'PO7', 'O1', 'Oz', 'O2', 'PO8', 'F3', 'Fz', 'F4', 'EOG'};
num_chans = length(chan_labels);
srate = 1000; 

% =================================================================
% 刺激検出パラメータ
ext_threshold = 2000; % 2200.037への立ち上がりを捉えるための閾値
stim_cooldown = 1000; % 検出後1000ms間は再検出を無視する (1.2秒周期に合わせる)
% =================================================================

for f = 1:length(files)
    file_path = fullfile(files(f).folder, files(f).name);
    [~, base_name, ~] = fileparts(files(f).name);
    fprintf('\n--- 処理中 (%d/%d): %s ---\n', f, length(files), base_name);
    
    % CSVの読み込み
    opts = detectImportOptions(file_path);
    opts.VariableNamesLine = 5; 
    opts.DataLines = [6, Inf]; 
    opts.PreserveVariableNames = true;
    
    try
        data_table = readtable(file_path, opts);
        data_table = rmmissing(data_table); 
        
        % 解析用EEGデータ（9ch）の抽出
        raw_data = table2array(data_table(:, chan_labels))'; 
        % 刺激検出用EXTデータの抽出
        ext_data = table2array(data_table(:, "EXT")); 
    catch ME
        fprintf('  警告: 読み込み失敗. %s\n', ME.message);
        continue;
    end
    
    % EEGLABデータセットの作成
    EEG = pop_importdata('dataformat', 'array', 'nbchan', num_chans, 'data', raw_data, ...
                         'srate', srate, 'pnts', size(raw_data, 2), 'xmin', 0);
    
    EEG.setname = base_name;
    for ch = 1:num_chans, EEG.chanlocs(ch).labels = chan_labels{ch}; end
    
    % NaNの物理除去
    bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
    if any(bad_samples)
        diff_nan = diff([0 bad_samples 0]);
        starts = find(diff_nan == 1);
        ends = find(diff_nan == -1) - 1;
        EEG = eeg_eegrej(EEG, [starts' ends']); 
    end
    
    % --- 【修正】クールダウン付き立ち上がり検出ロジック ---
    % 閾値を超えた全ての瞬間をまず抽出
    raw_stim_idx = find(diff(ext_data > ext_threshold) == 1) + 1; 
    
    stim_on_idx = [];
    if ~isempty(raw_stim_idx)
        last_idx = -stim_cooldown;
        for i = 1:length(raw_stim_idx)
            % 前回の検出地点から一定時間（1秒）以上経過している場合のみ採用
            if raw_stim_idx(i) - last_idx > stim_cooldown
                stim_on_idx = [stim_on_idx; raw_stim_idx(i)];
                last_idx = raw_stim_idx(i);
            end
        end
    end
    
    % --- 検出確認用グラフの自動保存 ---
    fig_check = figure('Name', ['Check Stimulus: ', base_name], 'Visible', 'off', 'Color', 'w');
    plot(ext_data); hold on;
    if ~isempty(stim_on_idx)
        plot(stim_on_idx, ext_data(stim_on_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end
    title(['Stimulus Onset Detection - ', base_name], 'Interpreter', 'none');
    ylabel('EXT Voltage');
    xlabel('Samples');
    legend('EXT Data', 'Detected Onset');
    grid on;
    
    save_path = fullfile(plot_save_dir, [base_name '_StimulusCheck.png']);
    exportgraphics(fig_check, save_path, 'Resolution', 150);
    close(fig_check);
    
    % イベントの登録
    EEG.event = [];
    if ~isempty(stim_on_idx)
        for i = 1:length(stim_on_idx)
            EEG.event(i).type = 'Saccade'; 
            EEG.event(i).latency = stim_on_idx(i);
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
        fprintf('  検出された刺激点数: %d\n', length(EEG.event));
    else
        fprintf('  警告: 刺激が検出されなかった.\n');
    end
    
    % 保存完了
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
    EEG = pop_saveset(EEG, 'filename', [base_name '.set'], 'filepath', output_dir, 'savemode', 'onefile');
end

fprintf('\n=== 修正版による変換と画像保存が完了した ===\n');