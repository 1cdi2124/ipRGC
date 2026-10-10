% ========================================================================
% ICA前後比較画像の単独実行スクリプト
% 使い方: MATLABコマンドウィンドウで「ICA_PLOT」と入力する。
% ========================================================================

input_dir = uigetdir('', 'サッケードOnset付き .set ファイルのフォルダを選択してください');
if input_dir == 0
    error('フォルダが選択されなかったため、処理を中断した。');
end

[parent_path, ~] = fileparts(input_dir);
main_output = fullfile(parent_path, 'ICA_PLOT');
iclabel_output = fullfile(main_output, 'ICLabel');
comparison_output = fullfile(main_output, 'ICA_before_after');
clean_set_output = fullfile(main_output, 'ICAcleaned_set');

output_dirs = {main_output, iclabel_output, comparison_output, clean_set_output};
for d = 1:numel(output_dirs)
    if ~exist(output_dirs{d}, 'dir')
        mkdir(output_dirs{d});
    end
end

file_list = dir(fullfile(input_dir, '*.set'));
if isempty(file_list)
    error('指定フォルダに .set ファイルが存在しない。');
end

% ICA_CWT_ERS_0702.m と揃えた設定
ica_chans = 1:8;
brain_threshold = 0.7;
plot_channel_labels = {'Oz'}; % 全チャネルなら {EEG.chanlocs(ica_chans).labels} に変更
event_type = 'Saccade';
window_ms = [-150 300];

if ~exist('eeglab', 'file')
    error('EEGLABのパスが通っていない。EEGLABを起動またはパス追加後に再実行する。');
end
if ~exist('pop_iclabel', 'file')
    error('ICLabelプラグインが見つからない。EEGLABのICLabelを導入してから再実行する。');
end

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
fprintf('\n=== ICA前後比較を開始: %d file(s) ===\n', numel(file_list));

for f = 1:numel(file_list)
    file_name = file_list(f).name;
    [~, base_name, ~] = fileparts(file_name);
    fprintf('\n--- 処理中 (%d/%d): %s ---\n', f, numel(file_list), file_name);

    try
        EEG = pop_loadset('filename', file_name, 'filepath', input_dir);

        bad_samples = any(isnan(EEG.data), 1) | any(isinf(EEG.data), 1);
        if any(bad_samples)
            diff_nan = diff([0 bad_samples 0]);
            starts = find(diff_nan == 1);
            ends = find(diff_nan == -1) - 1;
            EEG = eeg_eegrej(EEG, [starts' ends']);
            fprintf('  [Clean] NaN/Inf区間を削除した。\n');
        end

        if EEG.nbchan < max(ica_chans)
            error('ICA対象として必要な8 EEGチャネルがない（nbchan = %d）。', EEG.nbchan);
        end

        EEG.data = double(EEG.data);
        EEG = pop_chanedit(EEG, 'lookup', 'standard-10-5-cap385.elp');
        rng(3, 'twister');
        rng_state = rng;
        fprintf('  [RNG] Seed: %d\n', rng_state.Seed);

        evalc('EEG = pop_runica(EEG, ''icatype'', ''runica'', ''chanind'', ica_chans, ''extended'', 1);');
        EEG_for_label = pop_select(EEG, 'channel', ica_chans);
        evalc('EEG_for_label = pop_iclabel(EEG_for_label, ''default'');');

        classes = EEG_for_label.etc.ic_classification.ICLabel.classes;
        classifications = EEG_for_label.etc.ic_classification.ICLabel.classifications;
        brain_idx = find(strcmp(classes, 'Brain'), 1, 'first');
        if isempty(brain_idx)
            error('ICLabelの分類結果にBrainクラスが見つからない。');
        end
        auto_remove_ics = find(classifications(:, brain_idx) < brain_threshold);
        fprintf('  [ICLabel] 除去IC: %s\n', mat2str(auto_remove_ics));

        save_iclabel_report(EEG_for_label, classes, classifications, ...
            fullfile(iclabel_output, [base_name '_ICLabel.png']));

        EEG_before_ica = EEG;
        if ~isempty(auto_remove_ics)
            EEG = pop_subcomp(EEG, auto_remove_ics, 0);
        end
        EEG_after_ica = EEG;

        for c = 1:numel(plot_channel_labels)
            plot_channel = plot_channel_labels{c};
            try
                output_file = save_ica_before_after_figure( ...
                    EEG_before_ica, EEG_after_ica, comparison_output, base_name, ...
                    plot_channel, 'EOG', event_type, window_ms, auto_remove_ics);
                fprintf('  [Plot] %s\n', output_file);
            catch ME_plot
                warning('比較図の出力に失敗 (%s): %s', plot_channel, ME_plot.message);
            end
        end

        pop_saveset(EEG_after_ica, 'filename', [base_name '_ICAclean.set'], ...
            'filepath', clean_set_output, 'savemode', 'onefile');
        close all;

    catch ME
        warning('ファイル %s の処理をスキップ: %s', file_name, ME.message);
        close all;
    end
end

fprintf('\n=== ICA前後比較が完了した ===\n出力先: %s\n', main_output);


function save_iclabel_report(EEG_for_label, classes, classifications, output_file)
    num_ics = size(classifications, 1);
    fig = figure('Units', 'pixels', 'Position', [50 50 1600 1000], ...
        'Color', 'w', 'Visible', 'off');
    [spec, freqs] = spectopo(EEG_for_label.icaact, 0, EEG_for_label.srate, 'plot', 'off');

    for i = 1:num_ics
        subplot(4, 4, 2*i-1);
        topoplot(EEG_for_label.icawinv(:, i), EEG_for_label.chanlocs, ...
            'electrodes', 'on', 'style', 'both', 'headrad', 0.5);
        [max_prob, max_idx] = max(classifications(i, :));
        title(sprintf('IC %d: %s (%.1f%%)', i, classes{max_idx}, max_prob * 100), ...
            'FontSize', 16, 'FontWeight', 'bold');

        subplot(4, 4, 2*i);
        plot(freqs, spec(i, :), 'LineWidth', 2.0);
        xlim([1 100]); grid on;
        title(sprintf('IC %d Spectrum', i), 'FontSize', 14);
        set(gca, 'FontSize', 11, 'FontWeight', 'bold');
    end

    if exist('exportgraphics', 'file')
        exportgraphics(fig, output_file, 'Resolution', 200);
    else
        print(fig, output_file, '-dpng', '-r200');
    end
    close(fig);
end


function output_file = save_ica_before_after_figure(EEG_before, EEG_after, output_dir, ...
        base_name, channel_label, eog_label, event_type, window_ms, removed_ics)

    if ~ismatrix(EEG_before.data) || ~ismatrix(EEG_after.data)
        error('エポック化前の連続EEGを入力として想定している。');
    end
    if EEG_before.srate ~= EEG_after.srate || EEG_before.pnts ~= EEG_after.pnts
        error('ICA前後のEEGのサンプリングレートまたは長さが一致していない。');
    end

    eeg_idx = find_channel(EEG_before, channel_label);
    eog_idx = find_channel(EEG_before, eog_label);
    if isempty(eeg_idx)
        error('表示対象の電極 %s が見つからない。', channel_label);
    end

    event_types = arrayfun(@(e) event_type_to_char(e.type), EEG_before.event, 'UniformOutput', false);
    event_idx = find(strcmpi(event_types, event_type));
    if isempty(event_idx)
        error('イベント %s が見つからない。', event_type);
    end

    start_offset = round(window_ms(1) / 1000 * EEG_before.srate);
    end_offset = round(window_ms(2) / 1000 * EEG_before.srate);
    offsets = start_offset:end_offset;
    time_ms = offsets / EEG_before.srate * 1000;

    event_latencies = round([EEG_before.event(event_idx).latency]);
    is_valid = event_latencies + start_offset >= 1 & event_latencies + end_offset <= EEG_before.pnts;
    event_latencies = event_latencies(is_valid);
    if isempty(event_latencies)
        error('指定した範囲を切り出せるイベントがない。');
    end

    sample_matrix = bsxfun(@plus, event_latencies(:), offsets);
    before_trace = epoch_mean(EEG_before.data, eeg_idx, sample_matrix);
    after_trace = epoch_mean(EEG_after.data, eeg_idx, sample_matrix);
    removed_trace = before_trace - after_trace;
    if isempty(eog_idx)
        eog_trace = [];
    else
        eog_trace = epoch_mean(EEG_before.data, eog_idx, sample_matrix);
    end

    safe_base_name = regexprep(base_name, '[\\/:*?"<>|]', '_');
    output_file = fullfile(output_dir, ...
        sprintf('%s_ICA_before_after_%s.png', safe_base_name, channel_label));
    fig = figure('Color', 'w', 'Units', 'pixels', ...
        'Position', [80 80 1800 1200], 'Visible', 'off');

    subplot(3, 1, 1);
    if isempty(eog_trace)
        text(0.5, 0.5, 'EOG channel was not found', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 18);
        axis off;
    else
        plot(time_ms, eog_trace, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.8);
        grid on; xline(0, '--k', 'Saccade onset', 'LabelVerticalAlignment', 'bottom');
        xlim(window_ms); ylabel('EOG (µV)');
        title(sprintf('サッケード開始に整列した平均EOG（N = %d trials）', numel(event_latencies)), ...
            'FontSize', 18, 'FontWeight', 'bold');
        set(gca, 'FontSize', 14, 'LineWidth', 1.0);
    end

    common_limit = 1.1 * max(abs([before_trace, after_trace]));
    if common_limit == 0 || isnan(common_limit), common_limit = 1; end

    subplot(3, 1, 2);
    plot(time_ms, before_trace, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.8); hold on;
    plot(time_ms, after_trace, 'Color', [0.00 0.45 0.74], 'LineWidth', 2.0);
    grid on; xline(0, '--k', 'Saccade onset', 'LabelVerticalAlignment', 'bottom');
    xlim(window_ms); ylim([-common_limit common_limit]);
    ylabel(sprintf('%s (µV)', channel_label));
    title('ICA前後の平均波形（同一スケール）', 'FontSize', 18, 'FontWeight', 'bold');
    legend({'ICA前', 'ICA後'}, 'Location', 'best', 'FontSize', 14);
    set(gca, 'FontSize', 14, 'LineWidth', 1.0);

    subplot(3, 1, 3);
    plot(time_ms, removed_trace, 'Color', [0.85 0.33 0.10], 'LineWidth', 2.0);
    grid on; xline(0, '--k', 'Saccade onset', 'LabelVerticalAlignment', 'bottom');
    xlim(window_ms); ylabel(sprintf('%s (µV)', channel_label));
    xlabel('Time from saccade onset (ms)');
    title('ICAで除去された寄与（ICA前 − ICA後）', 'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'FontSize', 14, 'LineWidth', 1.0);

    if isempty(removed_ics)
        removed_ics_text = 'none';
    else
        removed_ics_text = strtrim(sprintf('%d ', removed_ics));
    end
    sgtitle(sprintf('%s | %s | Removed ICs: %s', ...
        strrep(base_name, '_', '\\_'), channel_label, removed_ics_text), ...
        'FontSize', 22, 'FontWeight', 'bold');

    if exist('exportgraphics', 'file')
        exportgraphics(fig, output_file, 'Resolution', 300);
    else
        print(fig, output_file, '-dpng', '-r300');
    end
    close(fig);
end


function chan_idx = find_channel(EEG, wanted_label)
    chan_idx = find(strcmpi({EEG.chanlocs.labels}, wanted_label), 1, 'first');
end


function event_type = event_type_to_char(raw_event_type)
    if ischar(raw_event_type)
        event_type = raw_event_type;
    elseif isstring(raw_event_type)
        event_type = char(raw_event_type);
    elseif isnumeric(raw_event_type)
        event_type = num2str(raw_event_type);
    else
        event_type = char(string(raw_event_type));
    end
end


function mean_trace = epoch_mean(data, chan_idx, sample_matrix)
    trials = reshape(double(data(chan_idx, sample_matrix(:))), size(sample_matrix));
    mean_trace = mean(trials, 1, 'omitnan');
end
