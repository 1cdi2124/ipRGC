% ================================================================= %
% lambda.m実行後，これで解析
% 目的：全条件に加え、各条件単独の3x3グラフを生成する。
% ================================================================= %

clear; close all;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- 1. ディレクトリ設定 ---
input_dir = uigetdir('', '解析対象の.setファイルが含まれるフォルダを選択してください');
if input_dir == 0, error('処理を中断しました。'); end

[parent_path, ~] = fileparts(input_dir);
main_output = fullfile(parent_path, '04_ERP');
if ~exist(main_output, 'dir'), mkdir(main_output); end

% --- 2. 解析パラメータの設定 ---
target_event = 'Saccade';        % サッケードマーカー
epoch_window = [-0.3 0.4];       % エポック区間
baseline_window = [-200 -100];   % ベースライン区間
plot_window = [-50 200];         % 描画区間

cond_keys = {'0O', '80O', '160O'}; 
cond_colors = {'k', 'b', 'r'}; 

file_list = dir(fullfile(input_dir, '*.set'));
if isempty(file_list), error('フォルダ内に.setファイルがありません。'); end

% --- 3. データの堅牢なスキャンとグループ化 ---
fprintf('ファイルをスキャンし、条件ごとのグループ化を行います...\n');
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
fprintf('検出された被験者数: %d\n', length(subjects));

% --- 4. データの読み込みと加算平均 ---
for sub_idx = 1:length(subjects)
    subj = subjects{sub_idx};
    fprintf('\n被験者 [%s] の処理中...\n', subj);
    
    conds_for_subj = fieldnames(group_struct.(subj));
    erp_storage = struct();
    time_points = [];
    chan_labels = {};
    num_channels = 0;
    
    for c = 1:length(conds_for_subj)
        c_id = conds_for_subj{c};
        f_name = group_struct.(subj).(c_id).f_name;
        orig_key = group_struct.(subj).(c_id).orig_key;
        
        fprintf('  -> 条件 %s ロード中: %s\n', orig_key, f_name);
        
        try
            EEG = pop_loadset('filename', f_name, 'filepath', input_dir);
            EEG_epoch = pop_epoch(EEG, {target_event}, epoch_window, 'newname', 'temp', 'epochinfo', 'yes');
            EEG_epoch = pop_rmbase(EEG_epoch, baseline_window);
            
            erp_storage.(c_id).data = mean(EEG_epoch.data, 3);
            erp_storage.(c_id).orig_key = orig_key;
            
            if isempty(time_points)
                time_points = EEG_epoch.times;
                num_channels = EEG_epoch.nbchan;
                chan_labels = {EEG_epoch.chanlocs.labels};
            end
        catch ME
            fprintf('    [エラー] %s の処理に失敗: %s\n', orig_key, ME.message);
        end
    end
    
    if isempty(fieldnames(erp_storage))
        fprintf('  警告: 有効なデータがロードできません。スキップします。\n');
        continue;
    end
    
    loaded_conds = fieldnames(erp_storage);
    
    % --- 5. 3x3サマリーグラフの描画（オーバーレイ ＋ 各条件単独） ---
    % 描画モードの定義（1つ目はオーバーレイ、以降はロードされた各条件）
    plot_modes = [{'Overlay'}; loaded_conds];
    
    chans_per_fig = 9; % 3x3
    num_figs = ceil(num_channels / chans_per_fig);
    
    for mode_idx = 1:length(plot_modes)
        current_mode = plot_modes{mode_idx};
        
        % 出力時のファイル名とタイトルの表示名を設定
        if strcmp(current_mode, 'Overlay')
            mode_display = 'Overlay';
        else
            mode_display = erp_storage.(current_mode).orig_key;
        end
        
        for fig_idx = 1:num_figs
            fig_summary = figure('Name', sprintf('%s - %s Summary Fig %d', subj, mode_display, fig_idx), ...
                         'Color', 'w', 'Position', [50, 50, 1600, 1000], 'Visible', 'off');
            
            start_ch = (fig_idx - 1) * chans_per_fig + 1;
            end_ch = min(fig_idx * chans_per_fig, num_channels);
            
            plot_idx = 1;
            for ch = start_ch:end_ch
                subplot(3, 3, plot_idx);
                legend_lines = [];
                legend_labels = {};
                
                % 描画する条件の決定（オーバーレイなら全て、それ以外なら単独）
                if strcmp(current_mode, 'Overlay')
                    conds_to_plot = loaded_conds;
                else
                    conds_to_plot = {current_mode};
                end
                
                for c = 1:length(conds_to_plot)
                    c_id = conds_to_plot{c};
                    orig_key = erp_storage.(c_id).orig_key;
                    
                    c_idx = find(strcmp(orig_key, cond_keys));
                    if ~isempty(c_idx), col = cond_colors{c_idx(1)}; else, col = 'g'; end
                    
                    h_plot = plot(time_points, erp_storage.(c_id).data(ch, :), 'LineWidth', 1.8, 'Color', col);
                    hold on;
                    legend_lines = [legend_lines, h_plot];
                    legend_labels = [legend_labels, {orig_key}];
                end
                
                xline(0, 'r--', 'Event', 'LineWidth', 1);
                yline(0, 'k-', 'LineWidth', 0.5);
                xlim(plot_window);
                xlabel('Time (ms)');
                ylabel('Amplitude (\muV)');
                title(sprintf('%s', chan_labels{ch}), 'Interpreter', 'none', 'FontSize', 12);
                grid on;
                
                if plot_idx == 1
                    legend(legend_lines, legend_labels, 'Location', 'best');
                end
                plot_idx = plot_idx + 1;
            end
            
            sgtitle(sprintf('%s | Mode: %s | Baseline [%d to %d ms]', subj, mode_display, baseline_window(1), baseline_window(2)), ...
                    'FontWeight', 'bold', 'FontSize', 16, 'Interpreter', 'none');
            
            save_name_summary = sprintf('%s_%s_ERP.png', subj, mode_display);
            exportgraphics(fig_summary, fullfile(main_output, save_name_summary), 'Resolution', 300);
            close(fig_summary);
        end
    end
    
    mat_save_path = fullfile(main_output, sprintf('%s_ERP_data.mat', subj));
    save(mat_save_path, 'erp_storage', 'time_points', 'chan_labels');
    fprintf('  -> 完了: 4パターンの画像を保存しました。\n');
end
fprintf('\n全処理完了。結果は %s に保存されました。\n', main_output);