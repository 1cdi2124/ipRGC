% =================================================================
% サッカード自動検出・バッチ処理スクリプト
% 処理フロー:
% 1. 生データの読み込みとクリーニング
% 2. 検出用の一時複製データ(EEG_filt)の作成とフィルタリング(0.5-30Hz)
% 3. 一時データを用いた6段階の品質ゲート評価とサッカード振幅頂点の割り出し
% 4. 割り出したイベント情報を元の生データ(EEG)へ登録
% 5. 生データ(.set)の保存 (CWTなどの広帯域解析へ接続可能)
% =================================================================

input_dir = uigetdir('', '生データ(CSV)が含まれるフォルダを選択してください');
if input_dir == 0
    error('フォルダが選択されなかったため, 処理を中断した.'); 
end

[parent_dir, ~, ~] = fileparts(input_dir);
output_dir = fullfile(parent_dir, '02_サッケード自動検出');
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end

files = dir(fullfile(input_dir, '*.csv'));
if isempty(files)
    error('指定されたディレクトリにCSVファイルが存在しない.'); 
end
fprintf('合計 %d 個のCSVファイルを処理する.\n', length(files));

if ~exist('eeglab', 'file')
    error('EEGLABのパスが通っていない. パスを追加してから再実行すること.');
end
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- 解析パラメータの定義 ---
chan_labels = {'PO7', 'O1', 'Oz', 'O2', 'PO8', 'F3', 'Fz', 'F4', 'EOG'};
num_chans = length(chan_labels);
eog_chan_idx = 9;
eeg_chan_indices = 1:8;
srate = 1000;

locutoff = 0.5;
hicutoff = 30.0;
threshold_z = 1.8; % 動的閾値
cooldown_ms = 250;
cooldown_samples = round((cooldown_ms / 1000) * srate);

% 単位補正済み品質ゲート閾値
MAX_EEG_P2P_THRESHOLD = 200.0; % 200 uV
MIN_STEP_AMP_THRESHOLD = 130.0; % 130 uV

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
        fprintf('  エラー: 読み込み失敗. スキップする. 詳細: %s\n', ME.message);
        continue;
    end
    
    % --- ステップ1: 生データをメモリに読み込む ---
    EEG = pop_importdata('dataformat', 'array', 'nbchan', num_chans, 'data', raw_data, ...
                         'srate', srate, 'pnts', size(raw_data, 2), 'xmin', 0);
    EEG.setname = base_name;
    for ch = 1:num_chans
        EEG.chanlocs(ch).labels = chan_labels{ch}; 
    end
    
    % 物理的な欠損値(NaN, Inf)のパージ (生データに対して適用)
    bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
    if any(bad_samples)
        diff_nan = diff([0 bad_samples 0]);
        starts = find(diff_nan == 1);
        ends = find(diff_nan == -1) - 1;
        EEG = eeg_eegrej(EEG, [starts' ends']); 
        fprintf('  [Clean] データ内の物理的な欠損区間を削除した.\n');
    end
    
    % --- ステップ2: 生データを複製し,複製したデータに対してのみフィルターをかける ---
    EEG_filt = EEG;
    
    try
        EEG_filt = pop_eegfiltnew(EEG_filt, 'locutoff', locutoff, 'hicutoff', hicutoff, 'channels', 1:num_chans);
        fprintf('  [Filter] 検出用一時データに %g-%g Hz のバンドパスフィルタを適用した.\n', locutoff, hicutoff);
    catch ME
        fprintf('  警告: フィルタリングに失敗. 生データで続行する. 詳細: %s\n', ME.message);
    end
    
    % --- ステップ3: フィルター適用済みのデータを使ってサッカードの位置を割り出す ---
    eog_data_filt = EEG_filt.data(eog_chan_idx, :); 
    eeg_data_filt = EEG_filt.data(eeg_chan_indices, :);
    
    eog_vel = [0, diff(eog_data_filt)] * srate; 
    abs_vel = abs(eog_vel);
    
    mean_vel = mean(abs_vel);
    std_vel = std(abs_vel);
    dyn_threshold = mean_vel + threshold_z * std_vel;
    fprintf('  [Threshold] 動的閾値: %.2f (Z=%.2f)\n', dyn_threshold, threshold_z);
    
    [candidate_peaks_vel, candidate_peaks_idx] = findpeaks(abs_vel, 'MinPeakHeight', dyn_threshold, 'MinPeakDistance', cooldown_samples);
    fprintf('  検出された候補数: %d\n', length(candidate_peaks_idx));
    
    valid_peaks = [];
    global_eog_std = std(eog_data_filt);
    
    if ~isempty(candidate_peaks_vel)
        q1 = prctile(candidate_peaks_vel, 25);
        q3 = prctile(candidate_peaks_vel, 75);
        iqr_val = q3 - q1;
        median_vel = median(candidate_peaks_vel);
        outlier_threshold = median_vel + 3 * iqr_val;
    else
        outlier_threshold = inf;
    end
    
    for i = 1:length(candidate_peaks_idx)
        peak = candidate_peaks_idx(i);
        
        if abs_vel(peak) > outlier_threshold
            continue;
        end
        
        win_mean = round(0.15 * srate); 
        gap = round(0.02 * srate);
        pre_start = max(1, peak - win_mean - gap);
        pre_end = max(1, peak - gap);
        post_start = min(length(eog_data_filt), peak + gap);
        post_end = min(length(eog_data_filt), peak + win_mean + gap);
        
        pre_mean = mean(eog_data_filt(pre_start:pre_end));
        post_mean = mean(eog_data_filt(post_start:post_end));
        step_amp = abs(post_mean - pre_mean);
        if step_amp < MIN_STEP_AMP_THRESHOLD
            continue;
        end
        
        win_overlap_start = round(0.1 * srate);
        win_overlap_end = round(0.3 * srate);
        start_eeg = max(1, peak - win_overlap_start);
        end_eeg = min(size(eeg_data_filt, 2), peak + win_overlap_end);
        
        ptp_per_chan = max(eeg_data_filt(:, start_eeg:end_eeg), [], 2) - min(eeg_data_filt(:, start_eeg:end_eeg), [], 2);
        max_eeg_ptp = max(ptp_per_chan);
        if max_eeg_ptp > MAX_EEG_P2P_THRESHOLD
            continue;
        end
        
        win_100ms = round(0.1 * srate);
        start_idx = max(1, peak - win_100ms);
        end_idx = min(length(abs_vel), peak + win_100ms);
        window_vel = abs_vel(start_idx:end_idx);
        primary_height = abs_vel(peak);
        
        [local_peaks_vel, local_peaks_locs] = findpeaks(window_vel, 'MinPeakProminence', 0.4 * primary_height);
        secondary_too_high = false;
        for lp = 1:length(local_peaks_locs)
            global_lp = start_idx + local_peaks_locs(lp) - 1;
            if abs(global_lp - peak) > 2
                secondary_too_high = true;
                break;
            end
        end
        if secondary_too_high
            continue;
        end
        
        win_200ms = round(0.2 * srate);
        end_idx_200 = min(length(eog_data_filt), peak + win_200ms);
        post_eog_window = eog_data_filt(peak:end_idx_200);
        if std(post_eog_window) >= 1.5 * global_eog_std
            continue;
        end
        
        win_300ms = round(0.3 * srate);
        end_idx_300 = min(length(eog_vel), peak + win_300ms);
        post_vel_window = eog_vel(peak:end_idx_300);
        
        kernel = ones(1, 25) / 25.0;
        smooth_vel = conv(post_vel_window, kernel, 'valid');
        deadband = 0.05 * primary_height;
        smooth_vel(abs(smooth_vel) < deadband) = 0;
        
        non_zero_vel = smooth_vel(smooth_vel ~= 0);
        if ~isempty(non_zero_vel)
            sign_changes = sum(diff(sign(non_zero_vel)) ~= 0);
        else
            sign_changes = 0;
        end
        
        if sign_changes > 4
            continue;
        end
        
        valid_peaks = [valid_peaks, peak];
    end
    
    fprintf('  有効なサッカード数: %d\n', length(valid_peaks));
    
    % --- ステップ4: イベント情報を「元の生データ(EEG)」に書き込む ---
    EEG.event = [];
    if ~isempty(valid_peaks)
        for i = 1:length(valid_peaks)
            p = valid_peaks(i);
            
            search_start = p;
            search_end = min(length(eog_data_filt), p + cooldown_samples);
            trigger_vel = eog_vel(p);
            
            if trigger_vel > 0
                [~, local_peak_rel_idx] = max(eog_data_filt(search_start:search_end));
            else
                [~, local_peak_rel_idx] = min(eog_data_filt(search_start:search_end));
            end
            
            amp_peak_idx = search_start + local_peak_rel_idx - 1;
            
            % オリジナルのEEG構造体にイベントを追加
            EEG.event(i).type = 'Saccade';
            EEG.event(i).latency = amp_peak_idx;
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end
    
    % 検証用グラフの作成 (プロットにはフィルタ適用済みのデータを使用)
    fig_check = figure('Name', ['Saccade Detection: ', base_name], 'Visible', 'off');
    
    subplot(2,1,1);
    plot(EEG_filt.times, abs_vel, 'Color', [0.2 0.4 0.6]); hold on;
    yline(dyn_threshold, 'r', 'LineWidth', 1.5, 'Label', sprintf('Threshold (Z=%.2f)', threshold_z));
    if ~isempty(valid_peaks)
        plot(EEG_filt.times(valid_peaks), abs_vel(valid_peaks), 'ro', 'MarkerFaceColor', 'r');
    end
    title(sprintf('EOG Velocity and Valid Peaks (Filtered Data) - %s', strrep(base_name, '_', '\_'))); 
    ylabel('Velocity (uV/s)'); grid on;
    
    subplot(2,1,2);
    plot(EEG_filt.times, eog_data_filt, 'Color', [0.1 0.1 0.1]); hold on;
    if ~isempty(EEG.event)
        for p = 1:length(EEG.event)
            line([EEG.times(EEG.event(p).latency) EEG.times(EEG.event(p).latency)], ylim, 'Color', 'r', 'LineWidth', 1.2);
        end
    end
    title('Filtered EOG Amplitude with Detected Saccade Amplitude Peaks');
    ylabel('Amplitude (uV)'); xlabel('Time (ms)'); grid on;
    
    save_fig_filename = fullfile(output_dir, [base_name, '.png']);
    saveas(fig_check, save_fig_filename);
    close(fig_check);
    
    % --- ステップ5: 生データ(イベントマーカー付き)を保存 ---
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0); 
    save_filename = [base_name '.set'];
    EEG = pop_saveset(EEG, 'filename', save_filename, 'filepath', output_dir, 'savemode', 'onefile');
end
fprintf('\n=== 全データの検出および【生データへのイベント登録】が完了した ===\n');
