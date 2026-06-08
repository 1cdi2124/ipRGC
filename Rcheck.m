% ================================================================= %
% サッケード持続時間 vs ERP振幅 相関解析プログラム（修正版）
% 目的：眼球運動の持続時間と脳波振幅の相関を計算し、
%       「遅延の正体」がアーティファクトか神経現象かを見極める。
% ================================================================= %

clear; close all;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- 1. 設定 ---
input_dir = uigetdir('', '解析済み.setファイルが含まれるフォルダを選択');
if input_dir == 0, error('処理を中断しました。'); end

parent_dir = fileparts(input_dir);
main_output = fullfile(parent_dir, '06_Correlation_Analysis');
if ~exist(main_output, 'dir'), mkdir(main_output); end

% 解析用パラメータ
target_chan = 'Oz'; % 解析対象チャネル
time_window = [0.05 0.25]; % Saccade_Start基準での期待区間（50-250ms）

% 【重要】最新の抽出アルゴリズムで登録したイベント名に合わせてください
start_marker = 'Saccade_Start';
end_marker   = 'Saccade_End';

file_list = dir(fullfile(input_dir, '*.set'));
group_struct = struct();

% --- 2. ファイルスキャン ---
for s = 1:length(file_list)
    f_name = file_list(s).name;
    [~, base_name, ~] = fileparts(f_name);
    
    % 被験者名の抽出（アンダースコアの最初をIDとする）
    parts = strsplit(base_name, '_');
    subj_name = matlab.lang.makeValidName(parts{1});
    
    if ~isfield(group_struct, subj_name)
        group_struct.(subj_name) = struct('files', {{}});
    end
    group_struct.(subj_name).files{end+1} = f_name;
end
subjects = fieldnames(group_struct);

% --- 3. 解析ループ ---
for sub_idx = 1:length(subjects)
    subj = subjects{sub_idx};
    fprintf('\n被験者 [%s] の相関解析を開始...\n', subj);
    
    all_durations = [];
    all_amplitudes = [];
    
    for f_idx = 1:length(group_struct.(subj).files)
        f_name = group_struct.(subj).files{f_idx};
        EEG = pop_loadset('filename', f_name, 'filepath', input_dir);
        
        % イベントリストから Saccade_Start と次の Saccade (End) を探す
        % エラー回避のためインデックスを慎重に操作
        events = {EEG.event.type};
        event_lats = [EEG.event.latency];
        
        for i = 1:length(events)-1
            if strcmp(events{i}, start_marker) && strcmp(events{i+1}, end_marker)
                start_lat = event_lats(i);
                end_lat = event_lats(i+1);
                
                % 持続時間の算出
                dur_ms = (end_lat - start_lat) / EEG.srate * 1000;
                
                % ノイズ除去（極端なサッケードを除外）
                if dur_ms > 20 && dur_ms < 200 
                    % ERP振幅の取得（Ozチャネル）
                    chan_idx = find(strcmp({EEG.chanlocs.labels}, target_chan));
                    if isempty(chan_idx), continue; end
                    
                    s_idx = round(start_lat + time_window(1)*EEG.srate);
                    e_idx = round(start_lat + time_window(2)*EEG.srate);
                    
                    % 境界チェック
                    s_idx = max(1, s_idx);
                    e_idx = min(EEG.pnts, e_idx);
                    
                    amp = max(abs(EEG.data(chan_idx, s_idx:e_idx))); % ピーク振幅
                    
                    all_durations(end+1) = dur_ms;
                    all_amplitudes(end+1) = amp;
                end
            end
        end
    end
    
    % --- 4. 相関解析とプロット ---
    if length(all_durations) > 10
        [R, P] = corrcoef(all_durations, all_amplitudes);
        
        fig = figure('Name', sprintf('%s Correlation', subj), 'Visible', 'off');
        scatter(all_durations, all_amplitudes, 'filled', 'MarkerFaceAlpha', 0.5);
        lsline; % 回帰直線を描画
        
        xlabel('Saccade Duration (ms)');
        ylabel(['ERP Amplitude (' target_chan ')']);
        title(sprintf('%s: R=%.3f, p=%.3e', subj, R(1,2), P(1,2)));
        grid on;
        
        save_path = fullfile(main_output, sprintf('%s_Correlation.png', subj));
        exportgraphics(fig, save_path, 'Resolution', 300);
        close(fig);
        
        fprintf('  被験者 %s: 相関係数 R=%.3f (p=%.3e)\n', subj, R(1,2), P(1,2));
    else
        fprintf('  被験者 %s: 十分な試行数がありませんでした。\n', subj);
    end
end
fprintf('\n解析完了。散布図は %s に保存されました。\n', main_output);