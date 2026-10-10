% ================================================================= %
% lambda.m実行後，これで解析
% 目的：全条件に加え、各条件単独の3x3グラフを生成する。
% さらに、全被験者の加算平均（Grand Average）も算出して出力する。
% ERP縦軸: 被験者ごとにEEG全ch・全条件で共通範囲を自動計算（EOGは個別自動）。
% Grand Averageも、平均波形からEEG共通範囲を別途計算する。
% ================================================================= %

clear; close all;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% --- 1. ディレクトリ設定 ---
% 必要な関数（pae_select_subject_input_dirs等）があるフォルダのパスを追加（どこからでも実行できるようにするため）
addpath('C:\Users\meizh\学校法人東海大学\高雄研 データ共有 - 光瞳孔反射研究プロジェクト\02瞳孔研究2025\05Phantom Array Effect\04 解析実行コード');

selected_root = uigetdir('', ...
    '02_サッケード自動検出、被験者フォルダ、または被験者一覧フォルダを選択');
if isequal(selected_root, 0), error('処理を中断しました。'); end
addpath(fileparts(mfilename('fullpath')), '-end');
input_dirs = pae_select_subject_input_dirs( ...
    selected_root, '02_サッケード自動検出', '.set');
fprintf('対象被験者フォルダ: %d 件\n', numel(input_dirs));

% --- 2. 解析パラメータの設定 ---
target_event = 'Saccade';        % サッケードマーカー
epoch_window = [-0.3 0.4];       % エポック区間
baseline_window = [-200 -100];   % ベースライン区間
plot_window = [-50 200];         % 描画区間
erp_axis_padding_fraction = 0.10; % EEG縦軸: 表示区間の最大絶対振幅に10%の余白
erp_eog_label_patterns = {'EOG'}; % 部分一致・大文字小文字不問。HEOG/VEOG/EOG1等を除外
% EEGは各被験者の全条件・全EEG chから0中心の対称範囲を決定（1 µV単位に切り上げ）。
% 同一被験者のOverlay・条件単独で共通。別被験者とGrand Averageでは別範囲。
% EOGは共通範囲へ含めず各グラフの自動縦軸を使う。ERP計算・保存値は変更しない。
validateattributes(erp_axis_padding_fraction, {'numeric'}, ...
    {'real','finite','scalar','nonnegative'}, mfilename, 'erp_axis_padding_fraction');
erp_eog_label_patterns = string(erp_eog_label_patterns(:));
if isempty(erp_eog_label_patterns) || any(ismissing(string(erp_eog_label_patterns))) || ...
        any(strlength(string(erp_eog_label_patterns)) == 0)
    error('lambdaERP:InvalidEOGPatterns','EOG判定パターンは空でない電極名文字列を指定してください。');
end

% 'NoGreen' と、それに対応する色（'g'）を追加
cond_keys = {'0O', '80O', '160O', 'NoGreen'}; 
cond_colors = {'k', 'b', 'r', 'g'}; 

% --- 全被験者平均（Grand Average）保存用変数 ---
grand_erp = struct();
grand_time_points = [];
grand_chan_labels = {};
grand_num_channels = 0;

for subject_index = 1:numel(input_dirs)
    input_dir = char(input_dirs(subject_index));
    [parent_path, ~] = fileparts(input_dir);
    main_output = fullfile(parent_path, '04_ERP');
    if ~exist(main_output, 'dir'), mkdir(main_output); end
    file_list = dir(fullfile(input_dir, '*.set'));
    fprintf('\n=== 被験者 %d/%d: %s（SET %d件）===\n', ...
        subject_index, numel(input_dirs), parent_path, numel(file_list));
    
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
            
            % 「緑なし」という文字列が含まれている場合は NoGreen に置換
            if contains(cond_str, '緑なし')
                cond_str = 'NoGreen';
            end
            
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
        [erp_ylim_uv,erp_yticks_uv,erp_eog_mask] = calculate_erp_plot_scale( ...
            erp_storage,time_points,chan_labels,plot_window, ...
            erp_eog_label_patterns,erp_axis_padding_fraction);
        
        % --- Grand Average 用のデータ蓄積 ---
        if isempty(grand_time_points)
            grand_time_points = time_points;
            grand_chan_labels = chan_labels;
            grand_num_channels = num_channels;
        end
        for c = 1:length(loaded_conds)
            c_id = loaded_conds{c};
            if ~isfield(grand_erp, c_id)
                grand_erp.(c_id).data_sum = erp_storage.(c_id).data;
                grand_erp.(c_id).count = 1;
                grand_erp.(c_id).orig_key = erp_storage.(c_id).orig_key;
            else
                grand_erp.(c_id).data_sum = grand_erp.(c_id).data_sum + erp_storage.(c_id).data;
                grand_erp.(c_id).count = grand_erp.(c_id).count + 1;
            end
        end
        
        % --- 5. 3x3サマリーグラフの描画（オーバーレイ ＋ 各条件単独） ---
        plot_modes = [{'Overlay'}; loaded_conds];
        chans_per_fig = 9; % 3x3
        num_figs = ceil(num_channels / chans_per_fig);
        
        for mode_idx = 1:length(plot_modes)
            current_mode = plot_modes{mode_idx};
            
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
                    if erp_eog_mask(ch)
                        ylim('auto'); yticks('auto');
                    else
                        ylim(erp_ylim_uv); yticks(erp_yticks_uv);
                    end
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
        save(mat_save_path, 'erp_storage', 'time_points', 'chan_labels', ...
            'erp_ylim_uv','erp_yticks_uv','erp_eog_mask', ...
            'erp_axis_padding_fraction','erp_eog_label_patterns');
        fprintf('  -> 完了: グラフ画像を保存しました。\n');
    end
end
fprintf('\n全 %d 被験者の処理が完了しました。\n', numel(input_dirs));

% ================================================================= %
% --- 6. 全被験者平均（Grand Average）の算出とグラフ出力 ---
% ================================================================= %
if ~isempty(fieldnames(grand_erp))
    fprintf('\n=== 全被験者平均 (Grand Average) の算出と出力 ===\n');
    ga_output = fullfile(selected_root, '04_ERP_GrandAverage');
    if ~exist(ga_output, 'dir'), mkdir(ga_output); end
    
    ga_storage = struct();
    ga_conds = fieldnames(grand_erp);
    
    % 平均の計算
    for c = 1:length(ga_conds)
        c_id = ga_conds{c};
        ga_storage.(c_id).data = grand_erp.(c_id).data_sum / grand_erp.(c_id).count;
        ga_storage.(c_id).orig_key = grand_erp.(c_id).orig_key;
        fprintf('  -> 条件 %s: %d 被験者のデータを平均\n', ga_storage.(c_id).orig_key, grand_erp.(c_id).count);
    end
    
    [erp_ylim_uv,erp_yticks_uv,erp_eog_mask] = calculate_erp_plot_scale( ...
        ga_storage,grand_time_points,grand_chan_labels,plot_window, ...
        erp_eog_label_patterns,erp_axis_padding_fraction);
    % 描画設定
    plot_modes_ga = [{'Overlay'}; ga_conds];
    chans_per_fig_ga = 9;
    num_figs_ga = ceil(grand_num_channels / chans_per_fig_ga);
    
    for mode_idx = 1:length(plot_modes_ga)
        current_mode = plot_modes_ga{mode_idx};
        
        if strcmp(current_mode, 'Overlay')
            mode_display = 'Overlay';
        else
            mode_display = ga_storage.(current_mode).orig_key;
        end
        
        for fig_idx = 1:num_figs_ga
            fig_summary_ga = figure('Name', sprintf('Grand Average - %s Summary Fig %d', mode_display, fig_idx), ...
                         'Color', 'w', 'Position', [50, 50, 1600, 1000], 'Visible', 'off');
            
            start_ch = (fig_idx - 1) * chans_per_fig_ga + 1;
            end_ch = min(fig_idx * chans_per_fig_ga, grand_num_channels);
            
            plot_idx = 1;
            for ch = start_ch:end_ch
                subplot(3, 3, plot_idx);
                legend_lines = [];
                legend_labels = {};
                
                if strcmp(current_mode, 'Overlay')
                    conds_to_plot = ga_conds;
                else
                    conds_to_plot = {current_mode};
                end
                
                for c = 1:length(conds_to_plot)
                    c_id = conds_to_plot{c};
                    orig_key = ga_storage.(c_id).orig_key;
                    
                    c_idx = find(strcmp(orig_key, cond_keys));
                    if ~isempty(c_idx), col = cond_colors{c_idx(1)}; else, col = 'g'; end
                    
                    h_plot = plot(grand_time_points, ga_storage.(c_id).data(ch, :), 'LineWidth', 2.0, 'Color', col);
                    hold on;
                    legend_lines = [legend_lines, h_plot];
                    legend_labels = [legend_labels, {orig_key}];
                end
                
                xline(0, 'r--', 'Event', 'LineWidth', 1);
                yline(0, 'k-', 'LineWidth', 0.5);
                xlim(plot_window);
                if erp_eog_mask(ch)
                    ylim('auto'); yticks('auto');
                else
                    ylim(erp_ylim_uv); yticks(erp_yticks_uv);
                end
                xlabel('Time (ms)');
                ylabel('Amplitude (\muV)');
                title(sprintf('%s', grand_chan_labels{ch}), 'Interpreter', 'none', 'FontSize', 12);
                grid on;
                
                if plot_idx == 1
                    legend(legend_lines, legend_labels, 'Location', 'best');
                end
                plot_idx = plot_idx + 1;
            end
            
            sgtitle(sprintf('Grand Average | Mode: %s | Baseline [%d to %d ms]', mode_display, baseline_window(1), baseline_window(2)), ...
                    'FontWeight', 'bold', 'FontSize', 16, 'Interpreter', 'none');
            
            save_name_summary_ga = sprintf('GrandAverage_%s_ERP.png', mode_display);
            exportgraphics(fig_summary_ga, fullfile(ga_output, save_name_summary_ga), 'Resolution', 300);
            close(fig_summary_ga);
        end
    end
    
    % Grand AverageのデータをMATファイルとして保存
    mat_save_path_ga = fullfile(ga_output, 'GrandAverage_ERP_data.mat');
    erp_storage = ga_storage; % 互換性のため変数名を寄せる
    time_points = grand_time_points;
    chan_labels = grand_chan_labels;
    save(mat_save_path_ga, 'erp_storage', 'time_points', 'chan_labels', ...
        'erp_ylim_uv','erp_yticks_uv','erp_eog_mask', ...
        'erp_axis_padding_fraction','erp_eog_label_patterns');
    fprintf('  -> Grand Average のグラフ画像とMATファイルを保存しました。\n');
else
    fprintf('\nGrand Average を計算するデータがありませんでした。\n');
end
fprintf('すべての処理が完了しました。\n');

function [limits,ticks,eog_mask] = calculate_erp_plot_scale(storage,times,labels,window,eog_patterns,padding)
% EEG全ch・全条件の表示区間だけから、その被験者/平均の共通範囲を求める。
    eog_mask = contains(upper(string(labels(:))),upper(string(eog_patterns(:))));
    limits = []; ticks = [];
    if all(eog_mask), return; end % EOGのみなら全グラフを自動調整する。
    use_time = times>=window(1) & times<=window(2);
    if ~any(use_time)
        error('lambdaERP:EmptyPlotWindow','表示区間にERPの時間点がありません。');
    end
    peak = 0; has_finite = false;
    conditions = fieldnames(storage);
    for index = 1:numel(conditions)
        data = storage.(conditions{index}).data;
        if ~isequal(size(data),[numel(labels),numel(times)])
            error('lambdaERP:ERPShape','条件 %s のERPと電極/時間のサイズが一致しません。',conditions{index});
        end
        values = double(data(~eog_mask,use_time));
        values = values(isfinite(values));
        if ~isempty(values)
            peak = max(peak,max(abs(values)));
            has_finite = true;
        end
    end
    if ~has_finite
        warning('lambdaERP:NoFiniteEEG','表示区間に有限なEEG値がないため、縦軸を仮に±1 µVにします。');
    end
    half_range = max(1,ceil(peak*(1+padding)));
    limits = [-half_range half_range];
    ticks = linspace(limits(1),limits(2),5);
end
