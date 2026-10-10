% =================================================================
% CWT解析・定量化一貫パイプライン（SPSS統計統合・シード追跡版）
% 目的: 高解像度CWT出力、シード値の動的取得とタイトル焼き付け、
%       8周波数帯×20ms刻みの厳密定量化、SPSS用一括統合マスターCSVの自動生成
% =================================================================

% --- ディレクトリ設定 ---
selected_root = uigetdir('', ...
    '03_手動サッケード、被験者フォルダ、または被験者一覧フォルダを選択');
if isequal(selected_root, 0), error('処理を中断した.'); end
addpath(fileparts(mfilename('fullpath')), '-end');
input_dirs = pae_select_subject_input_dirs( ...
    selected_root, '03_手動サッケード', '.set');
fprintf('対象被験者フォルダ: %d 件\n', numel(input_dirs));

% =================================================================
% --- 解析パラメータの設定 ---
% =================================================================
ica_chans = 1:8; 
epoch_window = [-0.3, 0.4]; 
baseline_window = [-250, -100]; 
freq_range = [0 200]; 
brain_threshold = 0.6;

% 高解像度（粒）の設定
n_freqs = 100;     % 周波数ステップ数（縦の細かさ）
n_timesout = 400;  % 時間ステップ数（横の細かさ）

% 色の強弱を調整
ersp_limit = 1.7; 
itc_limit = 0.55;   % ITCのコントラスト強調
% =================================================================

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultTextInterpreter', 'tex'); 
set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
set(groot, 'DefaultLegendInterpreter', 'tex');

for subject_index = 1:numel(input_dirs)
input_dir = char(input_dirs(subject_index));
[parent_path, ~] = fileparts(input_dir);
main_output = fullfile(parent_path, '04_CWT');
if ~exist(main_output, 'dir'), mkdir(main_output); end
plot_output = fullfile(main_output, 'ICLabel');
if ~exist(plot_output, 'dir'), mkdir(plot_output); end
file_list = dir(fullfile(input_dir, '*.set'));
fprintf('\n=== 被験者 %d/%d: %s（SET %d件）===\n', ...
    subject_index, numel(input_dirs), parent_path, numel(file_list));

% SPSS用マスターデータは被験者ごとに初期化する
spss_master_cell = {}; 

try
    for f = 1:length(file_list)
        file_name = file_list(f).name;
        [~, base_name, ~] = fileparts(file_name);
        fprintf('\n--- 処理開始 (%d/%d): %s ---\n', f, length(file_list), file_name);
        
        EEG = pop_loadset('filename', file_name, 'filepath', input_dir);
        
        % 1. 物理的異常区間除去
        bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
        if any(bad_samples)
            diff_nan = diff([0 bad_samples 0]);
            starts = find(diff_nan == 1);
            ends = find(diff_nan == -1) - 1;
            EEG = eeg_eegrej(EEG, [starts' ends']); 
            fprintf('  [Clean] 内部のNaN/Inf異常区間を物理的に切り詰めました.\n');
        end
        
        % 2. ICA ＆ ICLabel (倍精度化)
        EEG.data = double(EEG.data); 
        EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp');
        
        fprintf('  Removing 1 channel(s)...\n'); 
        
% =================================================================
        % 乱数シードの設定
        rng(308423991, 'twister'); % シード値
        % rng('shuffle');   % ランダム
        
        % 実際に適用されたシード値を取得して文字列に退避（再現性の確保）
        rng_state = rng;
        if isfield(rng_state, 'Seed')
            current_seed_str = num2str(rng_state.Seed);
        else
            current_seed_str = 'unknown';
        end
        fprintf('  [RNG] 使用中の乱数シード: %s\n', current_seed_str);
% =================================================================

        evalc('EEG = pop_runica(EEG, ''icatype'', ''runica'', ''chanind'', ica_chans, ''extended'', 1);');
        
        if exist('pop_iclabel', 'file')
            EEG_for_label = pop_select(EEG, 'channel', ica_chans);
            evalc('EEG_for_label = pop_iclabel(EEG_for_label, ''default'');');
            
            brain_idx = find(strcmp(EEG_for_label.etc.ic_classification.ICLabel.classes, 'Brain'));
            auto_remove_ics = find(EEG_for_label.etc.ic_classification.ICLabel.classifications(:, brain_idx) < brain_threshold);
            fprintf('  [AI] 脳活動判定: IC %s を除外対象に設定.\n', mat2str(auto_remove_ics));
            
            % --- ICLabel詳細レポート保存 ---
            try
                num_ics = length(ica_chans);
                fig_ic = figure('Units', 'pixels', 'Position', [50, 50, 1600, 1000], 'Color', 'w', 'WindowStyle', 'normal');
                clf(fig_ic);
                classes = EEG_for_label.etc.ic_classification.ICLabel.classes;
                classifications = EEG_for_label.etc.ic_classification.ICLabel.classifications;
                [spec, freqs_spec] = spectopo(EEG_for_label.icaact, 0, EEG_for_label.srate, 'plot', 'off');
                
                for i = 1:num_ics
                    subplot(4, 4, 2*i-1);
                    topoplot(EEG_for_label.icawinv(:, i), EEG_for_label.chanlocs, 'electrodes', 'on', 'style', 'both', 'headrad', 0.5);
                    [max_prob, max_idx] = max(classifications(i, :));
                    title(sprintf('IC %d: %s (%.1f%%)', i, classes{max_idx}, max_prob * 100), 'FontSize', 16, 'FontWeight', 'bold');
                    
                    subplot(4, 4, 2*i);
                    plot(freqs_spec, spec(i, :), 'LineWidth', 2.5);
                    xlim([1 100]); grid on;
                    title(['IC ' num2str(i) ' Spectrum'], 'FontSize', 14);
                    if i >= 7, xlabel('Freq (Hz)', 'FontSize', 12); end
                    set(gca, 'FontSize', 11, 'FontWeight', 'bold');
                end
                set(fig_ic, 'PaperPositionMode', 'auto');
                drawnow; pause(0.5); 
                save_path_ic = fullfile(plot_output, [base_name '_ICLabel.png']);
                if exist('exportgraphics', 'file')
                    exportgraphics(fig_ic, save_path_ic, 'Resolution', 150);
                else
                    saveas(fig_ic, save_path_ic);
                end
                close(fig_ic);
            catch
                fprintf('  警告: ICLabel画像の保存失敗.\n');
            end
            EEG.etc.ic_classification = EEG_for_label.etc.ic_classification;
        end
        
        % 3. 成分除去 ＆ エポック化
        if ~isempty(auto_remove_ics), EEG = pop_subcomp(EEG, auto_remove_ics, 0); end
        EEG = pop_epoch(EEG, {'Saccade'}, epoch_window, 'epochinfo', 'yes');
        EEG = pop_rmbase(EEG, baseline_window);
        
        % --- 全チャネル解析ループ ---
        all_target_chans = {EEG.chanlocs(1:8).labels};
        for c = 1:length(all_target_chans)
            this_chan = all_target_chans{c};
            chan_dir = fullfile(main_output, this_chan);
            if ~exist(chan_dir, 'dir'), mkdir(chan_dir); end
            
            chan_idx = find(strcmp({EEG.chanlocs.labels}, this_chan));
            
            % =================================================================
            % 1. ERSP & ITC の解析と保存 (上段・中段用)
            % =================================================================
            fig_ersp_itc = figure('Units', 'pixels', 'Position', [100, 100, 1200, 850], 'Visible', 'on', 'WindowStyle', 'normal');
            clf(fig_ersp_itc); 
            
            [ersp, itc, powbase, times, freqs] = newtimef(EEG.data(chan_idx, :, :), ...
                EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, 3, ...  
                'plotphase', 'off', 'padratio', 1, 'baseline', baseline_window, ...
                'freqs', freq_range, 'nfreqs', n_freqs, 'timesout', n_timesout, ...
                'erspmax', ersp_limit, 'itcmax', itc_limit, ...
                'plotersp', 'on', 'plotitc', 'on', 'title', ''); 
            
            % タイトル末尾にシード値を正確に焼き付ける
            title_str = [base_name ' : ' this_chan ' (ERSP / ITC / ERS) [Seed: ' current_seed_str ']'];
            sgtitle(title_str, 'FontSize', 17, 'FontWeight', 'bold', 'Interpreter', 'none');
            
            all_ax = findobj(fig_ersp_itc, 'Type', 'Axes');
            set(all_ax, 'FontSize', 12, 'FontWeight', 'bold'); 
            set(fig_ersp_itc, 'PaperPositionMode', 'auto');
            drawnow; pause(0.2); 
            
            temp_path_1 = fullfile(chan_dir, 'temp_ersp_itc.png');
            if exist('exportgraphics', 'file')
                exportgraphics(fig_ersp_itc, temp_path_1, 'Resolution', 150);
            else
                saveas(fig_ersp_itc, temp_path_1);
            end
            close(fig_ersp_itc); 
            
            % =================================================================
            % 2. ERS (絶対パワー) の解析と保存 (下段用)
            % =================================================================
            fig_ers = figure('Units', 'pixels', 'Position', [100, 100, 1200, 425], 'Visible', 'on', 'WindowStyle', 'normal');
            clf(fig_ers); 
            
            [ers, ~, ~, ~, ~] = newtimef(EEG.data(chan_idx, :, :), ...
                EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, 3, ...  
                'plotphase', 'off', 'padratio', 1, 'baseline', NaN, ...
                'freqs', freq_range, 'nfreqs', n_freqs, 'timesout', n_timesout, ...
                'plotersp', 'on', 'plotitc', 'off', 'title', ''); 
            
            % 1. カラーバー周辺の隠しラベルを置換
            all_cbs = findall(fig_ers, 'Type', 'colorbar');
            for i = 1:length(all_cbs)
                cb = all_cbs(i);
                if ischar(cb.YLabel.String) && (contains(cb.YLabel.String, 'ERSP') || contains(cb.YLabel.String, 'log10'))
                    cb.YLabel.String = 'ERS';
                    cb.YLabel.FontSize = 12;
                end
                if ischar(cb.Title.String) && (contains(cb.Title.String, 'ERSP') || contains(cb.Title.String, 'log10'))
                    cb.Title.String = 'ERS';
                end
            end
            
            % 2. 図内すべてのテキストオブジェクトを置換
            all_texts = findall(fig_ers, 'Type', 'text');
            for i = 1:length(all_texts)
                txt = all_texts(i).String;
                if ischar(txt) && (contains(txt, 'ERSP') || contains(txt, 'log10'))
                    all_texts(i).String = 'ERS(dB)';
                elseif iscell(txt)
                    for j = 1:length(txt)
                        if contains(txt{j}, 'ERSP') || contains(txt{j}, 'log10')
                            all_texts(i).String = 'ERS'; 
                            break;
                        end
                    end
                end
            end
            
            all_ax = findall(fig_ers, 'Type', 'Axes');
            set(all_ax, 'FontSize', 11, 'FontWeight', 'bold'); 
            set(fig_ers, 'PaperPositionMode', 'auto');
            drawnow; pause(0.2); 
            
            temp_path_2 = fullfile(chan_dir, 'temp_ers.png');
            if exist('exportgraphics', 'file')
                exportgraphics(fig_ers, temp_path_2, 'Resolution', 150);
            else
                saveas(fig_ers, temp_path_2);
            end
            close(fig_ers);
            
            % =================================================================
            % 3. 画像の結合 (3段構成の生成)
            % =================================================================
            try
                img1 = imread(temp_path_1);
                img2 = imread(temp_path_2);
                
                w1 = size(img1, 2);
                w2 = size(img2, 2);
                target_w = max(w1, w2);
                
                new_img1 = uint8(255 * ones(size(img1, 1), target_w, size(img1, 3)));
                new_img2 = uint8(255 * ones(size(img2, 1), target_w, size(img2, 3)));
                
                start1 = floor((target_w - w1) / 2) + 1;
                new_img1(:, start1 : start1+w1-1, :) = img1;
                
                start2 = floor((target_w - w2) / 2) + 1;
                new_img2(:, start2 : start2+w2-1, :) = img2;
                
                img_combined = [new_img1; new_img2];
                
                final_path = fullfile(chan_dir, [base_name '_' this_chan '.png']);
                imwrite(img_combined, final_path);
                
                delete(temp_path_1);
                delete(temp_path_2);
            catch ME
                fprintf('  警告: 画像の結合に失敗しました. ログ: %s\n', ME.message);
            end
            
            % =================================================================
            % 4. 厳密な8周波数帯域 × 20msブロック時間集約フェーズ (統計・SPSS用)
            % =================================================================
            % 20ms刻みの時間ブロック定義 (サッケード終了後 0ms から 200ms まで)
            time_edges = 0:20:200; 
            num_blocks = length(time_edges) - 1;
            
            % 8つの周波数帯域の定義（デルタからハイガンマ3分割）
            band_names = {'Delta', 'Theta', 'Alpha', 'Beta', 'LowGamma', 'HighGamma1', 'HighGamma2', 'HighGamma3'};
            band_ranges = [
                0.5,   3;   % Delta
                4,     7;   % Theta
                8,    13;   % Alpha
                14,   30;   % Beta
                31,   60;   % LowGamma
                61,  100;   % HighGamma (61~100)
                101, 150;   % HighGamma (101~150)
                151, 200    % HighGamma (151~200)
            ];
            num_bands = size(band_ranges, 1);
            
            % 個別保存用CSVの準備
            csv_headers = {'Frequency_Band'};
            for b = 1:num_blocks
                csv_headers{end+1} = sprintf('T_%d_to_%dms_ERSP', time_edges(b), time_edges(b+1));
                csv_headers{end+1} = sprintf('T_%d_to_%dms_ITC', time_edges(b), time_edges(b+1));
                csv_headers{end+1} = sprintf('T_%d_to_%dms_ERS', time_edges(b), time_edges(b+1));
            end
            
            individual_csv_cell = cell(num_bands + 1, length(csv_headers));
            individual_csv_cell(1, :) = csv_headers;
            
            % 定量結果を格納する構造体
            block_quant = struct();
            block_quant.meta.time_edges = time_edges;
            block_quant.meta.seed = current_seed_str;
            
            fprintf('  [Quant] 8周波数帯域 × 20ms集約を開始 (%s)...\n', this_chan);
            
            for b_idx = 1:num_bands
                b_name = band_names{b_idx};
                b_range = band_ranges(b_idx, :);
                
                % 該当する周波数インデックスの検索
                stat_f_idx = find(freqs >= b_range(1) & freqs <= b_range(2));
                if isempty(stat_f_idx)
                    [~, stat_f_idx] = min(abs(freqs - b_range(1))); 
                end
                
                % 時系列全体のトレース (周波数軸だけで平均化、時間400点の推移用)
                block_quant.(b_name).trace_ersp = mean(ersp(stat_f_idx, :), 1);
                block_quant.(b_name).trace_itc  = mean(abs(itc(stat_f_idx, :)), 1);
                block_quant.(b_name).trace_ers  = mean(ers(stat_f_idx, :), 1);
                
                individual_csv_cell{b_idx+1, 1} = b_name;
                col_ptr = 2;
                
                % SPSSマスター蓄積用行: [ファイル名, チャネル, 周波数帯域]
                spss_row = {base_name, this_chan, b_name};
                
                % 20msごとの時間ブロックループ
                for t_idx = 1:num_blocks
                    t_start = time_edges(t_idx);
                    t_end = time_edges(t_idx+1);
                    stat_t_idx = find(times >= t_start & times <= t_end);
                    
                    % 2次元ROI (周波数帯 × 時間ブロック) 内の平均値を算出
                    m_ersp = mean(mean(ersp(stat_f_idx, stat_t_idx), 2), 1);
                    m_itc  = mean(mean(abs(itc(stat_f_idx, stat_t_idx)), 2), 1);
                    m_ers  = mean(mean(ers(stat_f_idx, stat_t_idx), 2), 1);
                    
                    % 構造体へ格納
                    block_field = sprintf('block_%d_to_%dms', t_start, t_end);
                    block_quant.(b_name).(block_field).ersp = m_ersp;
                    block_quant.(b_name).(block_field).itc  = m_itc;
                    block_quant.(b_name).(block_field).ers  = m_ers;
                    
                    % 個別CSV用セルに格納
                    individual_csv_cell{b_idx+1, col_ptr}   = m_ersp;
                    individual_csv_cell{b_idx+1, col_ptr+1} = m_itc;
                    individual_csv_cell{b_idx+1, col_ptr+2} = m_ers;
                    col_ptr = col_ptr + 3;
                    
                    % SPSS一括用配列に格納
                    spss_row{end+1} = m_ersp;
                    spss_row{end+1} = m_itc;
                    spss_row{end+1} = m_ers;
                end
                
                % SPSS一括セル（全ファイル・チャネル・周波数帯の行を保持）に蓄積
                spss_master_cell(end+1, 1:length(spss_row)) = spss_row;
            end
            
            % 個別CSVをチャネルフォルダに書き出し
            csv_path = fullfile(chan_dir, [base_name '_' this_chan '_quant.csv']);
            writecell(individual_csv_cell, csv_path);
            
            % =================================================================
            % 5. データの保存 (MATデータ保存、block_quantを包含)
            % =================================================================
            save(fullfile(chan_dir, [base_name '_' this_chan '.mat']), 'ersp', 'itc', 'ers', 'powbase', 'times', 'freqs', 'block_quant');
            
            if c == 1
                pop_saveset(EEG, 'filename', [base_name '_Final.set'], 'filepath', chan_dir);
            end
        end
        close all;
    end
catch ME
    fprintf('  エラー発生: %s\n', ME.message);
end

% =================================================================
% SPSS専用：この被験者の全条件・全バンドをまとめたCSV
% =================================================================
if ~isempty(spss_master_cell)
    % ヘッダー行の生成
    time_edges = 0:20:200; 
    num_blocks = length(time_edges) - 1;
    master_headers = {'File_Name', 'Channel', 'Frequency_Band'};
    
    for b = 1:num_blocks
        master_headers{end+1} = sprintf('T_%d_%d_ERSP', time_edges(b), time_edges(b+1));
        master_headers{end+1} = sprintf('T_%d_%d_ITC', time_edges(b), time_edges(b+1));
        master_headers{end+1} = sprintf('T_%d_%d_ERS', time_edges(b), time_edges(b+1));
    end
    
    final_spss_table = [master_headers; spss_master_cell];
    master_csv_path = fullfile(main_output, 'SPSS.csv');
    writecell(final_spss_table, master_csv_path);
    fprintf('\n[Success] SPSS用CSVを出力しました:\n --> %s\n', master_csv_path);
else
    fprintf('\n[Warning] 処理データが蓄積されなかったため、CSVは生成されませんでした.\n');
end

fprintf('\n=== 全工程終了. 出力先: %s ===\n', main_output);
end

function s = num_idx_str(n)
    s = num2str(n);
end
