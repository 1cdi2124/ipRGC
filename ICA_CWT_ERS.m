% --- ディレクトリ設定 ---
input_dir = uigetdir('', '解析対象の.setファイルが含まれるフォルダを選択してください');
if input_dir == 0, error('処理を中断した.'); end

[parent_path, ~] = fileparts(input_dir); 
[~, parent_folder_name] = fileparts(parent_path); 

clean_name = strrep(parent_folder_name, '_', '');
%output_folder_name = [clean_name, '_CWT'];
output_folder_name = '04_CWT';
main_output = fullfile(parent_path, output_folder_name);

if ~exist(main_output, 'dir'), mkdir(main_output); end
plot_output = fullfile(main_output, 'ICLabel');
if ~exist(plot_output, 'dir'), mkdir(plot_output); end

file_list = dir(fullfile(input_dir, '*.set'));

% =================================================================
% --- 解析パラメータの設定 ---調整
ica_chans = 1:8; 
epoch_window = [-0.3, 0.4]; 
baseline_window = [-250, -100]; 
freq_range = [0 200]; 
brain_threshold = 0.5;

% 高解像度（粒）の設定
n_freqs = 100;     % 周波数ステップ数（縦の細かさ）
n_timesout = 400;  % 時間ステップ数（横の細かさ）

% 色の強弱を調整
%ersp_limit = 2.75;  % ERSPのコントラスト強調
ersp_limit = 1.7; 
itc_limit = 0.55;   % ITCのコントラスト強調
% =================================================================

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultTextInterpreter', 'tex'); 
set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
set(groot, 'DefaultLegendInterpreter', 'tex');

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
        rng(3, 'twister'); % シード値
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
            
            sgtitle([base_name ' : ' this_chan ' (ERSP / ITC / ERS)'], 'FontSize', 17, 'FontWeight', 'bold', 'Interpreter', 'none');
            
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
            % 高さを約半分にし、1段分のサイズに最適化
            fig_ers = figure('Units', 'pixels', 'Position', [100, 100, 1200, 425], 'Visible', 'on', 'WindowStyle', 'normal');
            clf(fig_ers); 
            
            % baselineをNaNに設定することで基線補正を回避し、対数絶対パワー(ERS)を算出
            [ers, ~, ~, ~, ~] = newtimef(EEG.data(chan_idx, :, :), ...
                EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, 3, ...  
                'plotphase', 'off', 'padratio', 1, 'baseline', NaN, ...
                'freqs', freq_range, 'nfreqs', n_freqs, 'timesout', n_timesout, ...
                'plotersp', 'on', 'plotitc', 'off', 'title', ''); 

            % 1. カラーバー周辺の隠しラベルを検索・置換
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

            % 2. 図内すべてのテキストオブジェクト（複数行セル配列含む）を検索・置換
            all_texts = findall(fig_ers, 'Type', 'text');
            for i = 1:length(all_texts)
                txt = all_texts(i).String;
                if ischar(txt) && (contains(txt, 'ERSP') || contains(txt, 'log10'))
                    all_texts(i).String = 'ERS(dB)';
                elseif iscell(txt)
                    % テキストが複数行（ERSP \n (10*log10...)）で構成されている場合
                    for j = 1:length(txt)
                        if contains(txt{j}, 'ERSP') || contains(txt{j}, 'log10')
                            all_texts(i).String = 'ERS'; % 複数行を1行の 'ERS' に統合して上書き
                            break;
                        end
                    end
                end
            end
            % -------------------------------------------------------------
            
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
                
                % 白背景のキャンバスを作成して中央揃えで結合
                new_img1 = uint8(255 * ones(size(img1, 1), target_w, size(img1, 3)));
                new_img2 = uint8(255 * ones(size(img2, 1), target_w, size(img2, 3)));
                
                start1 = floor((target_w - w1) / 2) + 1;
                new_img1(:, start1 : start1+w1-1, :) = img1;
                
                start2 = floor((target_w - w2) / 2) + 1;
                new_img2(:, start2 : start2+w2-1, :) = img2;
                
                img_combined = [new_img1; new_img2];
                
                final_path = fullfile(chan_dir, [base_name '_' this_chan '.png']);
                imwrite(img_combined, final_path);
                
                % 一時ファイルのクリーンアップ
                delete(temp_path_1);
                delete(temp_path_2);
            catch ME
                fprintf('  警告: 画像の結合に失敗しました. ログ: %s\n', ME.message);
            end

            % =================================================================
            % データの保存 (ERSP, ITC, ERS をすべて格納)
            % =================================================================
            save(fullfile(chan_dir, [base_name '_' this_chan '.mat']), 'ersp', 'itc', 'ers', 'powbase', 'times', 'freqs');
            
            if c == 1
                pop_saveset(EEG, 'filename', [base_name '_Final.set'], 'filepath', chan_dir);
            end
        end
        close all;
    end
catch ME
    fprintf('  エラー発生: %s\n', ME.message);
end

fprintf('\n=== 全工程終了. 出力先: %s ===\n', main_output);

function s = num_idx_str(n)
    s = num2str(n);
end