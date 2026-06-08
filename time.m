% ================================================================= %
% サッケード持続時間（行動データ）算出・比較コード
% mix.m実行後，これで開始と終了までの時間を計測
% 目的：0Hz, 80Hz, 160Hz条件におけるサッケード持続時間の平均値を計算し、
%       ERPで見られた「遅延」の原因が眼球運動の差であるかを検証する。
% ================================================================= %

clear; close all;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- 1. ディレクトリ設定 ---
input_dir = uigetdir('', '解析対象の.setファイルが含まれるフォルダを選択してください');
if input_dir == 0, error('処理を中断しました。'); end

[parent_path, ~] = fileparts(input_dir);
output_dir = fullfile(parent_path, '05_time');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% --- 2. 解析パラメータの設定 ---
% 【警告】必ずあなたのデータセット内の実際のイベント名に書き換えてください！
start_marker = 'Saccade_Start'; % サッケード開始のイベント名
end_marker   = 'Saccade_End';       % サッケード終了のイベント名

cond_keys = {'0O', '80O', '160O'};
cond_colors = {[0 0 0], [0 0 1], [1 0 0]}; % 0O=黒, 80O=青, 160O=赤

file_list = dir(fullfile(input_dir, '*.set'));
if isempty(file_list), error('フォルダ内に.setファイルがありません。'); end

% --- 3. データのグループ化 ---
group_struct = struct();
for s = 1:length(file_list)
    f_name = file_list(s).name;
    [~, base_name, ~] = fileparts(f_name);
    
    parts = strsplit(base_name, '_');
    if length(parts) >= 2
        subj_name = parts{1}; 
        cond_str  = parts{2}; 
        
        if ismember(cond_str, cond_keys)
            s_id = matlab.lang.makeValidName(subj_name);
            c_id = matlab.lang.makeValidName(cond_str);
            group_struct.(s_id).(c_id).f_name = f_name;
            group_struct.(s_id).(c_id).orig_key = cond_str;
        end
    end
end

subjects = fieldnames(group_struct);

% --- 4. 持続時間の計算 ---
for sub_idx = 1:length(subjects)
    subj = subjects{sub_idx};
    fprintf('\n====================================================\n');
    fprintf('被験者 [%s] のサッケード持続時間解析を開始\n', subj);
    fprintf('====================================================\n');
    
    conds_for_subj = fieldnames(group_struct.(subj));
    duration_stats = struct();
    
    for c = 1:length(conds_for_subj)
        c_id = conds_for_subj{c};
        f_name = group_struct.(subj).(c_id).f_name;
        orig_key = group_struct.(subj).(c_id).orig_key;
        
        try
            % データのロード（波形ではなくイベント情報のみを使用）
            EEG = pop_loadset('filename', f_name, 'filepath', input_dir);
            
            durations = [];
            % イベントリストから Start と End のペアを探す
            for i = 1:length(EEG.event)-1
                if strcmp(EEG.event(i).type, start_marker)
                    % 開始マーカーが見つかったら、直後の終了マーカーを探す
                    for j = i+1:min(i+5, length(EEG.event)) 
                        if strcmp(EEG.event(j).type, end_marker)
                            % レイテンシ（データポイント）の差からミリ秒を計算
                            dur_ms = (EEG.event(j).latency - EEG.event(i).latency) / EEG.srate * 1000;
                            
                            % 極端な異常値（例: 5ms以下、200ms以上など）はノイズとして除外
                            if dur_ms > 5 && dur_ms < 200
                                durations(end+1) = dur_ms;
                            end
                            break; % 次の Start を探すために内側のループを抜ける
                        end
                    end
                end
            end
            
            if isempty(durations)
                fprintf('  [警告] %s: 指定された Start (%s) と End (%s) のペアが見つかりませんでした。\n', orig_key, start_marker, end_marker);
                continue;
            end
            
            % 統計値の計算
            mean_dur = mean(durations);
            std_dur = std(durations);
            trial_num = length(durations);
            
            duration_stats.(c_id).mean = mean_dur;
            duration_stats.(c_id).std = std_dur;
            duration_stats.(c_id).orig_key = orig_key;
            duration_stats.(c_id).trial_num = trial_num;
            
            fprintf('  条件 %-5s: 平均 %5.2f ms (SD: %5.2f) | 有効試行数: %d\n', orig_key, mean_dur, std_dur, trial_num);
            
        catch ME
            fprintf('  [エラー] %s の処理に失敗: %s\n', orig_key, ME.message);
        end
    end
    
    % --- 5. 比較グラフの描画 ---
    loaded_conds = fieldnames(duration_stats);
    if length(loaded_conds) > 1
        fig = figure('Name', sprintf('%s - Saccade Durations', subj), 'Color', 'w', 'Position', [100, 100, 600, 500]);
        hold on;
        
        means = [];
        stds = [];
        labels = {};
        colors_to_use = [];
        
        for c = 1:length(cond_keys)
            c_id = matlab.lang.makeValidName(cond_keys{c});
            if isfield(duration_stats, c_id)
                means(end+1) = duration_stats.(c_id).mean;
                stds(end+1)  = duration_stats.(c_id).std;
                labels{end+1} = duration_stats.(c_id).orig_key;
                colors_to_use = [colors_to_use; cond_colors{c}];
            end
        end
        
        % 棒グラフの描画
        for i = 1:length(means)
            bar(i, means(i), 'FaceColor', colors_to_use(i,:), 'EdgeColor', 'k', 'LineWidth', 1.5, 'BarWidth', 0.6);
        end
        % エラーバー（標準偏差）の描画
        errorbar(1:length(means), means, stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
        
        set(gca, 'XTick', 1:length(means), 'XTickLabel', labels, 'FontSize', 14);
        ylabel('Saccade Duration (ms)', 'FontSize', 14);
        title(sprintf('Subject: %s | Mean Saccade Duration', subj), 'FontSize', 16, 'Interpreter', 'none');
        grid on;
        
        % 画像保存
        save_name = sprintf('%s_Saccade_Duration.png', subj);
        exportgraphics(fig, fullfile(output_dir, save_name), 'Resolution', 300);
        
        % 統計データの保存
        mat_save_path = fullfile(output_dir, sprintf('%s_Duration_Stats.mat', subj));
        save(mat_save_path, 'duration_stats');
    end
end

fprintf('\nすべての持続時間解析が完了しました。結果は %s に保存されています。\n', output_dir);