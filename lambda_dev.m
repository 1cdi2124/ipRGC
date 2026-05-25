% =================================================================
% サッケード有り(O)を解析する
% 頂点検出モデル: 検出点からcooldown期間内の絶対振幅最大値をイベントとする
% グラフ保存機能付き
% =================================================================

% 1. フォルダの選択
input_dir = uigetdir('', '生データ(CSV)が含まれるフォルダを選択してください');
if input_dir == 0, error('フォルダが選択されなかったため, 処理を中断した.'); end

% 保存先フォルダの作成
[parent_dir, ~, ~] = fileparts(input_dir);
output_dir = fullfile(parent_dir, '02_サッケード自動検出');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% 重複を避けるため, *.csv で一括取得
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
threshold = 40; % 速度の閾値
cooldown = 400; % 探索窓兼クールダウン(ms)
% =================================================================

for f = 1:length(files)
    file_path = fullfile(files(f).folder, files(f).name);
    [~, base_name, ~] = fileparts(files(f).name);
    fprintf('\n--- 処理中 (%d/%d): %s ---\n', f, length(files), base_name);
    
    opts = detectImportOptions(file_path);
    opts.VariableNamesLine = 5; 
    opts.DataLines = [6, Inf];  
    
    try
        data_table = readtable(file_path, opts);
        data_table = rmmissing(data_table); 
        raw_data = table2array(data_table(:, chan_labels))'; 
    catch ME
        fprintf('  警告: 読み込み失敗. %s\n', ME.message);
        continue;
    end
    
    EEG = pop_importdata('dataformat', 'array', 'nbchan', num_chans, 'data', raw_data, ...
                         'srate', srate, 'pnts', size(raw_data, 2), 'xmin', 0);
    EEG.setname = base_name;
    for ch = 1:num_chans, EEG.chanlocs(ch).labels = chan_labels{ch}; end
    
    bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
    if any(bad_samples)
        diff_nan = diff([0 bad_samples 0]);
        starts = find(diff_nan == 1);
        ends = find(diff_nan == -1) - 1;
        EEG = eeg_eegrej(EEG, [starts' ends']); 
        fprintf('  [Clean] 配列内の物理的なNaN区間を削除しました.\n');
    end

    % --- サッカード検出および頂点補正ロジック ---
    eog_data = EEG.data(9, :); 
    eog_diff = [0, diff(eog_data)]; 
    
    % 速度ベースでの初期検出
    initial_idx = find(abs(eog_diff) > threshold);
    
    peak_saccades = [];
    last_idx = -cooldown;
    
    for i = 1:length(initial_idx)
        if initial_idx(i) - last_idx > cooldown
            % 検出点からcooldown期間内(全サンプル)で振幅の絶対値が最大となる地点を探す
            search_start = initial_idx(i);
            search_end = min(initial_idx(i) + cooldown, length(eog_data));
            
            [~, local_peak_rel_idx] = max(abs(eog_data(search_start:search_end)));
            peak_idx = search_start + local_peak_rel_idx - 1;
            
            peak_saccades = [peak_saccades; peak_idx];
            last_idx = initial_idx(i); % 次の検出判定は元の検出点基準でクールダウン
        end
    end
    
    % --- グラフ作成と保存 ---
    fig_check = figure('Name', ['Check Peak & Threshold: ', base_name], 'Visible', 'off'); % 画面に表示せずバックグラウンドで描画
    
    % 上段: 速度(diff)と閾値の確認
    subplot(2,1,1);
    plot(abs(eog_diff)); hold on;
    yline(threshold, 'r', 'LineWidth', 1.5, 'Label', 'Current Threshold');
    title(['Velocity (|diff|) - ', base_name]); 
    ylabel('|diff|'); grid on;
    
    % 下段: 元の振幅と検出された頂点の確認
    subplot(2,1,2);
    plot(eog_data); hold on;
    if ~isempty(peak_saccades)
        for p = 1:length(peak_saccades)
            line([peak_saccades(p) peak_saccades(p)], ylim, 'Color', 'r', 'LineWidth', 1.2);
        end
    end
    title('EOG Amplitude with Peak Events (Red Lines)');
    ylabel('Amplitude'); xlabel('Samples'); grid on;
    
    % グラフを画像として保存
    save_fig_filename = fullfile(output_dir, [base_name, '_threshold_check.png']);
    saveas(fig_check, save_fig_filename);
    fprintf('  グラフを保存した: %s\n', save_fig_filename);
    
    % イベントの登録
    EEG.event = [];
    if ~isempty(peak_saccades)
        for i = 1:length(peak_saccades)
            EEG.event(i).type = 'Saccade';
            EEG.event(i).latency = peak_saccades(i);
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
        fprintf('  検出されたサッカード数: %d (頂点に配置完了)\n', length(EEG.event));
    else
        fprintf('  警告: サッカードが検出されなかった.\n');
    end
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
    save_filename = [base_name '.set'];
    EEG = pop_saveset(EEG, 'filename', save_filename, 'filepath', output_dir, 'savemode', 'onefile');
    
    close(fig_check);
end

fprintf('\n=== 頂点基準でのサンプリング, グラフ保存, および保存が完了した ===\n');