% =================================================================
% 生データに開始と終了どっちものイベント設置
% サッカード自動検出・バッチ処理スクリプト (Start & End 同時抽出・堅牢版)
% 信号処理(バンドパスフィルタ, 微分スケーリング, Zスコア動的閾値)
% 検出基準: 
%   1. Start: 速度(微分値)が動的閾値を超えた最初の瞬間
%   2. End  : Start後クールダウン期間内における運動方向に応じた振幅の頂点
% =================================================================

% 1. フォルダの選択
input_dir = uigetdir('', '生データ(CSV)が含まれるフォルダを選択してください');
if input_dir == 0
    error('フォルダが選択されなかったため, 処理を中断した.'); 
end

% 保存先フォルダの作成
[parent_dir, ~, ~] = fileparts(input_dir);
output_dir = fullfile(parent_dir, '02_サッケード自動検出');
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end

% CSVファイルの一括取得
files = dir(fullfile(input_dir, '*.csv'));
if isempty(files)
    error('指定されたディレクトリにCSVファイルが存在しない.'); 
end
fprintf('合計 %d 個のCSVファイルを処理する.\n', length(files));

% EEGLABの初期化とパス確認
if ~exist('eeglab', 'file')
    error('EEGLABのパスが通っていない. パスを追加してから再実行すること.');
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% =================================================================
% 解析パラメータの定義
% =================================================================
chan_labels = {'PO7', 'O1', 'Oz', 'O2', 'PO8', 'F3', 'Fz', 'F4', 'EOG'};
num_chans = length(chan_labels);
eog_chan_idx = 9; 
srate = 1000; 

% フィルタリングパラメータ
locutoff = 0.5; 
hicutoff = 30.0; 

% 検出アルゴリズムパラメータ
threshold_z = 1.85; 
cooldown_ms = 250; 
cooldown_samples = round((cooldown_ms / 1000) * srate);

% =================================================================
for f = 1:length(files)
    file_path = fullfile(files(f).folder, files(f).name);
    [~, base_name, ~] = fileparts(files(f).name);
    fprintf('\n--- 処理中 (%d/%d): %s ---\n', f, length(files), base_name);
    
    % CSVの読み込み処理
    opts = detectImportOptions(file_path);
    opts.VariableNamesLine = 5; 
    opts.DataLines = [6, Inf];  
    
    try
        data_table = readtable(file_path, opts);
        data_table = rmmissing(data_table); 
        raw_data = table2array(data_table(:, chan_labels))'; 
    catch ME
        fprintf('  エラー: 読み込み失敗. スキップする. 詳細: %s\n', ME.message);
        continue;
    end
    
    % EEGLABデータ構造への格納
    EEG = pop_importdata('dataformat', 'array', 'nbchan', num_chans, 'data', raw_data, ...
                         'srate', srate, 'pnts', size(raw_data, 2), 'xmin', 0);
    EEG.setname = base_name;
    for ch = 1:num_chans
        EEG.chanlocs(ch).labels = chan_labels{ch}; 
    end
    
    % 物理的な欠損値(NaN, Inf)のパージ
    bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
    if any(bad_samples)
        diff_nan = diff([0 bad_samples 0]);
        starts = find(diff_nan == 1);
        ends = find(diff_nan == -1) - 1;
        EEG = eeg_eegrej(EEG, [starts' ends']); 
        fprintf('  [Clean] データ内の物理的な欠損区間を削除した.\n');
    end
    
    % 帯域通過フィルタの適用
    try
        EEG = pop_eegfiltnew(EEG, 'locutoff', locutoff, 'hicutoff', hicutoff, 'channels', eog_chan_idx);
        fprintf('  [Filter] EOGに %g-%g Hz のバンドパスフィルタを適用.\n', locutoff, hicutoff);
    catch ME
        fprintf('  警告: フィルタリングに失敗. 生データで続行する. 詳細: %s\n', ME.message);
    end
    
    % 速度(微分)の計算とスケーリング
    eog_data = EEG.data(eog_chan_idx, :); 
    eog_vel = [0, diff(eog_data)] * srate; 
    abs_vel = abs(eog_vel);
    
    % Zスコアに基づく動的閾値の算出
    mean_vel = mean(abs_vel);
    std_vel = std(abs_vel);
    dyn_threshold = mean_vel + threshold_z * std_vel;
    fprintf('  [Threshold] 動的閾値: %.2f (Z=%.2f)\n', dyn_threshold, threshold_z);
    
    % --- Start & End 同時抽出ロジック ---
    initial_idx = find(abs_vel > dyn_threshold);
    
    % ペアを格納する配列
    start_saccades = [];
    end_saccades = [];
    
    last_idx = -cooldown_samples;
    
    for i = 1:length(initial_idx)
        if initial_idx(i) - last_idx > cooldown_samples
            % 1. Start地点の確定 (閾値を超えた瞬間)
            search_start = initial_idx(i);
            search_end = min(initial_idx(i) + cooldown_samples, length(eog_data));
            
            % 2. 閾値を超えた瞬間の速度の正負でサッカードの方向を判定
            trigger_vel = eog_vel(search_start);
            
            if trigger_vel > 0
                % 陽性方向への移動: 振幅の最大値を探索
                [~, local_peak_rel_idx] = max(eog_data(search_start:search_end));
            else
                % 陰性方向への移動: 振幅の最小値を探索
                [~, local_peak_rel_idx] = min(eog_data(search_start:search_end));
            end
            
            % 3. End地点の確定 (振幅の頂点)
            peak_idx = search_start + local_peak_rel_idx - 1;
            
            % 4. 1対1のペアとして配列に格納
            start_saccades = [start_saccades; search_start];
            end_saccades = [end_saccades; peak_idx];
            
            % 次の探索は今回のEnd点を基準にクールダウンを適用
            last_idx = peak_idx; 
        end
    end
    
    % --- グラフ作成と保存 (検証用視覚化) ---
    fig_check = figure('Name', ['Saccade Detection: ', base_name], 'Visible', 'off');
    
    % 上段: 速度と動的閾値 (Startを青丸でプロット)
    subplot(2,1,1);
    plot(EEG.times, abs_vel, 'Color', [0.2 0.4 0.6]); hold on;
    yline(dyn_threshold, 'r', 'LineWidth', 1.5, 'Label', sprintf('Threshold (Z=%.2f)', threshold_z));
    if ~isempty(start_saccades)
        plot(EEG.times(start_saccades), abs_vel(start_saccades), 'bo', 'MarkerFaceColor', 'b');
    end
    title(sprintf('EOG Velocity with Start Points (Blue) - %s', strrep(base_name, '_', '\_'))); 
    ylabel('Velocity (uV/s)'); grid on;
    
    % 下段: フィルタリング後のEOG振幅 (Startを青線、Endを赤線でプロット)
    subplot(2,1,2);
    plot(EEG.times, eog_data, 'Color', [0.1 0.1 0.1]); hold on;
    if ~isempty(end_saccades)
        for p = 1:length(end_saccades)
            line([EEG.times(start_saccades(p)) EEG.times(start_saccades(p))], ylim, 'Color', 'b', 'LineWidth', 1.0, 'LineStyle', '--');
            line([EEG.times(end_saccades(p)) EEG.times(end_saccades(p))], ylim, 'Color', 'r', 'LineWidth', 1.5);
        end
    end
    title('EOG Amplitude with Start (Blue dash) and End (Red solid)');
    ylabel('Amplitude (uV)'); xlabel('Time (ms)'); grid on;
    
    % 画像の書き出し
    save_fig_filename = fullfile(output_dir, [base_name, '.png']);
    saveas(fig_check, save_fig_filename);
    close(fig_check);
    
    % --- EEGLABへのイベント登録 ---
    EEG.event = [];
    event_idx = 1;
    
    if ~isempty(end_saccades) && ~isempty(start_saccades)
        % 安全装置: 数が一致しない場合はエラーで止める
        if length(start_saccades) ~= length(end_saccades)
            error('StartとEndの検出数が一致していません。アルゴリズムが破綻しています。');
        end
        
        for i = 1:length(end_saccades)
            % Startイベントの登録
            EEG.event(event_idx).type = 'Saccade_Start';
            EEG.event(event_idx).latency = start_saccades(i);
            event_idx = event_idx + 1;
            
            % Endイベントの登録
            EEG.event(event_idx).type = 'Saccade_End';
            EEG.event(event_idx).latency = end_saccades(i);
            event_idx = event_idx + 1;
        end
        
        EEG = eeg_checkset(EEG, 'eventconsistency');
        fprintf('  [Result] 検出完了. Start/Endペア数: %d組 (総イベント数: %d)\n', length(end_saccades), length(EEG.event));
    else
        fprintf('  [Result] 警告: サッカードは検出されなかった.\n');
    end
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
    save_filename = [base_name '.set'];
    EEG = pop_saveset(EEG, 'filename', save_filename, 'filepath', output_dir, 'savemode', 'onefile');
end
fprintf('\n=== 全データのサッカード Start & End 同時検出が完了した ===\n');