function master_manifest = run_set_cwt_band_timecourse(input_root)
%RUN_SET_CWT_BAND_TIMECOURSE Onset SETからCWT/ERSP帯域時系列を一括再計算する。
%　緑なし差分とバンド帯両方
%
% 実行例:
%   run_set_cwt_band_timecourse
%   run_set_cwt_band_timecourse("C:\...\04_模擬実験_80_160_CWT_1000Hz")
%   run_set_cwt_band_timecourse("C:\...\takahira_0427\03_手動サッケード")
%
% 入力ルート以下の「03_手動サッケード」フォルダを再帰検索し、全SETを処理する。
% 出力は各被験者フォルダの「06_バンド帯\各チャンネル」以下に保存する。
%
% 重要:
%   - 0 Hzは振動ではなくDCなので解析しない。最低周波数は0.5 Hz。
%   - 条件名の0Oは緑LED常灯の非PAE対照、緑なしは無灯の
%     眼球運動対照であり、別条件として扱う。
%   - 緑なしがあれば「各条件-緑なし」を追加出力する。
%   - Delta/Thetaは連続データ上でCWTを行ってからサッケード時刻に揃える。
%   - サッケード間隔が約0.66秒のデータでは、Delta/Thetaは隣接する
%     サッケードを含む探索的指標であり、単一サッケード反応とはみなさない。
%   - Alpha以上は長めのエポックでnewtimefを行い、-200～300 msを表示する。

    config = default_config();

    if nargin < 1 || strlength(string(input_root)) == 0
        selected = uigetdir('', ...
            ['実験ルート、被験者フォルダ、または' ...
             '03_手動サッケードを選択してください']);
        if isequal(selected, 0)
            error('run_set_cwt_band_timecourse:Cancelled', ...
                '処理を中断しました。');
        end
        input_root = string(selected);
    else
        input_root = string(input_root);
    end

    if ~isfolder(input_root)
        error('run_set_cwt_band_timecourse:MissingInputFolder', ...
            '入力フォルダが見つかりません: %s', input_root);
    end

    initialize_eeglab(config);
    manual_folders = find_manual_saccade_folders(input_root);
    if isempty(manual_folders)
        error('run_set_cwt_band_timecourse:NoManualSaccadeFolder', ...
            '入力ルート以下に03_手動サッケードが見つかりません: %s', ...
            input_root);
    end

    master_manifest = table();
    total_difference_outputs = 0;
    fprintf('\n=== SETからのCWT/ERSP帯域解析を開始 ===\n');
    fprintf('検索ルート: %s\n', input_root);
    fprintf('対象03_手動サッケード数: %d\n', numel(manual_folders));
    fprintf('表示範囲: %.0f～%.0f ms\n', ...
        config.display_window_ms(1), config.display_window_ms(2));

    for folder_index = 1:numel(manual_folders)
        manual_folder = manual_folders(folder_index);
        subject_folder = string(fileparts(manual_folder));
        output_root = fullfile(subject_folder, config.output_relative_path);
        if ~isfolder(output_root)
            mkdir(output_root);
        end

        set_list = dir(fullfile(manual_folder, '*.set'));
        fprintf('\n--- 入力フォルダ: %s ---\n', manual_folder);
        fprintf('SET数: %d / 出力: %s\n', numel(set_list), output_root);

        folder_results = cell(0, 1);
        for file_index = 1:numel(set_list)
            set_path = string(fullfile(set_list(file_index).folder, ...
                set_list(file_index).name));
            fprintf('\n[%d/%d] %s\n', file_index, numel(set_list), ...
                set_list(file_index).name);

            try
                result = process_one_set(set_path, output_root, config);
                folder_results{end + 1, 1} = result; %#ok<AGROW>
                row = make_manifest_row(manual_folder, set_path, ...
                    output_root, "完了", "", result);
                fprintf('  完了: unique Saccade=%d, ICA除去=%s\n', ...
                    result.UniqueSaccadeEvents, result.RemovedICs);
            catch ME
                warning('run_set_cwt_band_timecourse:FileFailed', ...
                    '%s の処理に失敗しました: %s', set_path, ME.message);
                row = make_manifest_row(manual_folder, set_path, ...
                    output_root, "失敗", string(ME.message), struct());
            end

            if isempty(master_manifest)
                master_manifest = row;
            else
                master_manifest = [master_manifest; row]; %#ok<AGROW>
            end
        end

        difference_manifest = create_green_difference_outputs( ...
            folder_results, output_root, config);
        if ~isempty(difference_manifest)
            writetable(difference_manifest, ...
                fullfile(output_root, 'SET_CWT_緑なし差分_manifest.csv'), ...
                'Encoding', 'UTF-8');
            total_difference_outputs = total_difference_outputs + ...
                sum(difference_manifest.Status == "完了");
        end

        folder_mask = master_manifest.ManualSaccadeFolder == manual_folder;
        writetable(master_manifest(folder_mask, :), ...
            fullfile(output_root, 'SET_CWT_manifest.csv'), ...
            'Encoding', 'UTF-8');
    end

    fprintf('\n=== SETからのCWT/ERSP帯域解析が終了しました ===\n');
    fprintf('完了SET: %d / 失敗SET: %d\n', ...
        sum(master_manifest.Status == "完了"), ...
        sum(master_manifest.Status == "失敗"));
    fprintf('緑なし差分ERSP出力: %dチャンネル\n', ...
        total_difference_outputs);
    fprintf(['注意: Delta/Thetaは隣接サッケードが重なるため、' ...
        '探索的結果として解釈してください。\n']);
end

function config = default_config()
    config.output_relative_path = "06_バンド帯_test";
    config.target_channels = ["F3", "F4", "Fz", "O1", ...
        "O2", "Oz", "PO7", "PO8"];
    config.event_type = "Saccade";
    config.event_duplicate_tolerance_samples = 1;

    config.display_window_ms = [-200 300];
    config.baseline_window_ms = [-300 -200];
    config.output_step_ms = 5;
    config.common_times_ms = ...
        config.display_window_ms(1):config.output_step_ms:...
        config.display_window_ms(2);

    config.low_frequencies_hz = 0.5:0.5:7;
    config.high_frequencies_hz = 8:1:200;
    config.wavelet_cycles = 3;
    config.high_epoch_window_s = [-1 1];
    config.high_timesout = 401;
    config.padratio = 1;

    config.band_names = [ ...
        "Delta"; "Theta"; "Alpha"; "Beta"; "LowGamma"; ...
        "HighGamma_61_100"; "HighGamma_101_150"; ...
        "HighGamma_151_200"];
    config.band_labels = [ ...
        "デルタ"; "シータ"; "アルファ"; "ベータ"; ...
        "ローガンマ"; "ハイガンマ1"; "ハイガンマ2"; ...
        "ハイガンマ3"];
    config.band_ranges = [ ...
          0.5,   3; 4,   7; 8,  13; 14,  30; ...
         31,    60; 61, 100; 101, 150; 151, 200];
    config.band_interpretation = [ ...
        "探索的（隣接サッケード重畳）"; ...
        "探索的（隣接サッケード重畳）"; ...
        "主解析"; "主解析"; "主解析"; "主解析"; ...
        "主解析"; "主解析"];

    % 商用電源周波数と高調波はCWT自体には残すが、帯域平均から除外する。
    config.exclude_line_noise_bins = true;
    config.line_frequencies_hz = [50 100 150 200];
    config.line_exclusion_half_width_hz = 1;

    % Brain<閾値ではなく、明確なアーチファクト確率だけで除去する。
    config.run_ica = true;
    config.ica_seed = 3;
    config.artifact_class_names = [ ...
        "Muscle", "Eye", "Heart", "Line Noise", "Channel Noise"];
    config.artifact_thresholds = [0.80 0.80 0.90 0.90 0.90];

    config.eeglab_fallback = ...
        "C:\Users\tatsuya\Downloads\Apps\eeglab2025.1.0";
end

function initialize_eeglab(config)
    if exist('eeglab', 'file') ~= 2 && isfolder(config.eeglab_fallback)
        addpath(config.eeglab_fallback);
    end
    if exist('eeglab', 'file') ~= 2
        error('run_set_cwt_band_timecourse:EEGLABNotFound', ...
            ['EEGLABがMATLABパスにありません。EEGLABを起動してから' ...
             '再実行してください。']);
    end
    if exist('newtimef', 'file') ~= 2 || exist('timefreq', 'file') ~= 2
        evalc('eeglab nogui');
    elseif exist('pop_iclabel', 'file') ~= 2
        % プラグインパスを含めて初期化する。
        evalc('eeglab nogui');
    end
end

function manual_folders = find_manual_saccade_folders(input_root)
    input_root = string(input_root);
    [~, selected_name] = fileparts(input_root);
    if startsWith(string(selected_name), "03") && ...
            contains(string(selected_name), "手動サッケード")
        manual_folders = input_root;
        return;
    end

    sets = dir(fullfile(input_root, '**', '*.set'));
    if isempty(sets)
        manual_folders = strings(0, 1);
        return;
    end
    folders = unique(string({sets.folder}), 'stable');
    keep = false(size(folders));
    for index = 1:numel(folders)
        [~, folder_name] = fileparts(folders(index));
        keep(index) = startsWith(string(folder_name), "03") && ...
            contains(string(folder_name), "手動サッケード");
    end
    manual_folders = folders(keep);
end

function result = process_one_set(set_path, output_root, config)
    [~, base_name] = fileparts(set_path);
    safe_base = sanitize_filename(base_name);
    qc_folder = fullfile(output_root, "ICA_QC");
    if ~isfolder(qc_folder)
        mkdir(qc_folder);
    end

    EEG = load_set_with_fdt_repair(set_path);
    if EEG.trials ~= 1
        error('run_set_cwt_band_timecourse:NotContinuous', ...
            '%s は連続データではありません（trials=%d）。', ...
            set_path, EEG.trials);
    end
    if abs(EEG.srate - 1000) > 1e-6
        warning('run_set_cwt_band_timecourse:UnexpectedSamplingRate', ...
            '%s のサンプリング周波数は%.3g Hzです。', set_path, EEG.srate);
    end

    [EEG, raw_saccades, duplicates_removed] = ...
        deduplicate_saccade_events(EEG, config.event_type, ...
            config.event_duplicate_tolerance_samples);

    bad_samples = any(~isfinite(EEG.data), 1);
    if any(bad_samples)
        changes = diff([false, bad_samples, false]);
        starts = find(changes == 1);
        stops = find(changes == -1) - 1;
        EEG = eeg_eegrej(EEG, [starts(:), stops(:)]);
        fprintf('  NaN/Inf区間を%d区間除去しました。\n', numel(starts));
    end

    target_indices = find_target_channels(EEG, config.target_channels);
    [EEG, ica_report] = run_ica_cleaning(EEG, target_indices, ...
        qc_folder, safe_base, config);

    event_latencies = get_event_latencies(EEG, config.event_type);
    if isempty(event_latencies)
        error('run_set_cwt_band_timecourse:NoSaccadeEvents', ...
            '%s にSaccadeイベントがありません。', set_path);
    end

    EEG_high = pop_epoch(EEG, {char(config.event_type)}, ...
        config.high_epoch_window_s, 'epochinfo', 'yes');
    if EEG_high.trials < 2
        error('run_set_cwt_band_timecourse:TooFewEpochs', ...
            '%s の有効エポックが不足しています。', set_path);
    end

    all_summary = table();
    all_timecourses = table();
    low_valid_event_count = NaN;
    cwt_output_files = strings(numel(config.target_channels), 1);
    band_output_files = strings(numel(config.target_channels), 1);

    for channel_index = 1:numel(config.target_channels)
        channel_name = config.target_channels(channel_index);
        continuous_channel_index = target_indices(channel_index);
        epoched_channel_index = find(strcmpi( ...
            string({EEG_high.chanlocs.labels}), channel_name), 1);

        fprintf('  CWT: %s ...\n', channel_name);
        [low_ersp, low_freqs, low_valid] = ...
            compute_low_frequency_ersp(EEG, continuous_channel_index, ...
                event_latencies, config);
        [high_ersp, high_freqs] = ...
            compute_high_frequency_ersp(EEG_high, ...
                epoched_channel_index, config);
        if isnan(low_valid_event_count)
            low_valid_event_count = low_valid;
        end

        freqs = [low_freqs(:); high_freqs(:)];
        times = config.common_times_ms;
        ersp = [low_ersp; high_ersp];

        [band_timecourses, band_summary, timecourse_table] = ...
            calculate_band_timecourses(freqs, ersp, times, config);

        channel_output = fullfile(output_root, channel_name);
        if ~isfolder(channel_output)
            mkdir(channel_output);
        end
        cwt_stem = safe_base + "_" + channel_name + "_CWT再計算";
        band_stem = safe_base + "_" + channel_name + ...
            "_SET再計算_band_timecourse";
        output_cwt_mat = fullfile(channel_output, cwt_stem + ".mat");
        output_cwt_png = fullfile(channel_output, cwt_stem + ".png");
        output_band_mat = fullfile(channel_output, band_stem + ".mat");
        output_band_png = fullfile(channel_output, band_stem + ".png");
        output_timecourse_csv = fullfile(channel_output, ...
            band_stem + ".csv");
        output_summary_csv = fullfile(channel_output, ...
            band_stem + "_summary.csv");
        cwt_output_files(channel_index) = output_cwt_mat;
        band_output_files(channel_index) = output_band_mat;

        source_set = set_path;
        analysis_config = config;
        save(output_cwt_mat, 'source_set', 'analysis_config', ...
            'ersp', 'freqs', 'times', 'low_ersp', 'low_freqs', ...
            'high_ersp', 'high_freqs', 'low_valid', '-v7');
        save(output_band_mat, 'source_set', 'analysis_config', ...
            'times', 'band_timecourses', 'band_summary', '-v7');
        writetable(timecourse_table, output_timecourse_csv, ...
            'Encoding', 'UTF-8');
        writetable(band_summary, output_summary_csv, ...
            'Encoding', 'UTF-8');
        create_band_figure(times, band_timecourses, band_summary, ...
            base_name, channel_name, output_band_png, config);
        create_cwt_figure(times, freqs, ersp, base_name, ...
            channel_name, output_cwt_png, config);

        summary_with_source = addvars(band_summary, ...
            repmat(channel_name, height(band_summary), 1), ...
            repmat(set_path, height(band_summary), 1), ...
            'Before', 1, ...
            'NewVariableNames', {'Channel', 'SourceSET'});
        timecourse_with_source = addvars(timecourse_table, ...
            repmat(channel_name, height(timecourse_table), 1), ...
            repmat(set_path, height(timecourse_table), 1), ...
            'Before', 1, ...
            'NewVariableNames', {'Channel', 'SourceSET'});

        if isempty(all_summary)
            all_summary = summary_with_source;
            all_timecourses = timecourse_with_source;
        else
            all_summary = [all_summary; summary_with_source]; %#ok<AGROW>
            all_timecourses = [all_timecourses; ...
                timecourse_with_source]; %#ok<AGROW>
        end
    end

    writetable(all_summary, ...
        fullfile(output_root, safe_base + "_全チャンネル帯域集計.csv"), ...
        'Encoding', 'UTF-8');
    writetable(all_timecourses, ...
        fullfile(output_root, safe_base + "_全チャンネル時系列.csv"), ...
        'Encoding', 'UTF-8');

    result.RawSaccadeEvents = raw_saccades;
    result.UniqueSaccadeEvents = numel(event_latencies);
    result.DuplicatesRemoved = duplicates_removed;
    result.HighEpochCount = EEG_high.trials;
    result.LowValidEventCount = low_valid_event_count;
    result.RemovedICs = string(mat2str(ica_report.RemovedICs));
    result.SourceSET = set_path;
    result.BaseName = string(base_name);
    result.ChannelNames = config.target_channels(:);
    result.CWTFiles = cwt_output_files;
    result.BandFiles = band_output_files;
end

function EEG = load_set_with_fdt_repair(set_path)
    set_path = string(set_path);
    [set_folder, set_base, set_extension] = fileparts(set_path);
    loaded = load(char(set_path), '-mat');
    if isfield(loaded, 'EEG')
        EEG = loaded.EEG;
    else
        EEG = loaded;
    end

    EEG.filepath = char(set_folder);
    EEG.filename = char(set_base + set_extension);
    if ~isnumeric(EEG.data)
        expected_fdt = fullfile(set_folder, set_base + ".fdt");
        stored_fdt = fullfile(set_folder, string(EEG.data));
        if isfile(expected_fdt)
            EEG.data = char(set_base + ".fdt");
            EEG.datfile = EEG.data;
        elseif isfile(stored_fdt)
            [~, stored_name, stored_extension] = fileparts(stored_fdt);
            EEG.data = char(stored_name + stored_extension);
            EEG.datfile = EEG.data;
        else
            candidates = dir(fullfile(set_folder, '*.fdt'));
            if isscalar(candidates)
                EEG.data = candidates(1).name;
                EEG.datfile = candidates(1).name;
            else
                error('run_set_cwt_band_timecourse:FDTNotFound', ...
                    '%s に対応するFDTを一意に決められません。', set_path);
            end
        end
        EEG.data = eeg_getdatact(EEG);
    end
    EEG.data = double(EEG.data);
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

function [EEG, raw_count, removed_count] = ...
        deduplicate_saccade_events(EEG, event_type, tolerance_samples)
    types = event_types_as_strings(EEG.event);
    target_indices = find(strcmpi(types, event_type));
    raw_count = numel(target_indices);
    if raw_count < 2
        removed_count = 0;
        return;
    end

    latencies = double([EEG.event(target_indices).latency]);
    [sorted_latencies, order] = sort(latencies);
    duplicate_sorted = [false, ...
        diff(sorted_latencies) <= tolerance_samples];
    remove_indices = target_indices(order(duplicate_sorted));
    EEG.event(remove_indices) = [];
    EEG = eeg_checkset(EEG, 'eventconsistency');
    removed_count = numel(remove_indices);
    if removed_count > 0
        fprintf('  重複Saccadeイベントを%d個除去しました。\n', removed_count);
    end
end

function target_indices = find_target_channels(EEG, target_channels)
    labels = string({EEG.chanlocs.labels});
    target_indices = zeros(1, numel(target_channels));
    for index = 1:numel(target_channels)
        found = find(strcmpi(labels, target_channels(index)), 1);
        if isempty(found)
            error('run_set_cwt_band_timecourse:MissingChannel', ...
                'チャンネル%sがSETにありません。', target_channels(index));
        end
        target_indices(index) = found;
    end
end

function [EEG, report] = run_ica_cleaning(EEG, target_indices, ...
        qc_folder, safe_base, config)
    report.RemovedICs = [];
    if ~config.run_ica
        return;
    end

    if exist('pop_runica', 'file') ~= 2
        warning('run_set_cwt_band_timecourse:RunICANotFound', ...
            'pop_runicaがないためICAを省略します。');
        return;
    end

    try
        if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
            EEG = pop_chanedit(EEG, ...
                'lookup', 'standard-10-5-cap385.elp');
        end
    catch ME
        warning('run_set_cwt_band_timecourse:ChannelLookupFailed', ...
            '標準座標の補完に失敗しました: %s', ME.message);
    end

    rng(config.ica_seed, 'twister');
    rank_value = rank(double(EEG.data(target_indices, :))');
    runica_options = {'icatype', 'runica', 'chanind', target_indices, ...
        'extended', 1};
    if rank_value < numel(target_indices)
        runica_options = [runica_options, {'pca', rank_value}];
    end
    fprintf('  ICA実行: EEG %d ch, rank=%d\n', ...
        numel(target_indices), rank_value);
    EEG = pop_runica(EEG, runica_options{:});

    if exist('pop_iclabel', 'file') ~= 2
        warning('run_set_cwt_band_timecourse:ICLabelNotFound', ...
            'ICLabelがないため成分自動除去を行いません。');
        return;
    end

    % ICLabelはEEG.icachansindを参照できるため、EOGを含む元のEEGを
    % pop_selectで縮小しない。縮小するとICA activation/weight mismatchの
    % 警告が生じる場合がある。
    EEG_for_label = EEG;
    evalc('EEG_for_label = pop_iclabel(EEG_for_label, ''default'');');
    classification = ...
        EEG_for_label.etc.ic_classification.ICLabel.classifications;
    classes = string( ...
        EEG_for_label.etc.ic_classification.ICLabel.classes);
    remove_mask = false(size(classification, 1), 1);
    removal_reason = strings(size(classification, 1), 1);

    for class_index = 1:numel(config.artifact_class_names)
        class_name = config.artifact_class_names(class_index);
        probability_column = find(strcmpi(classes, class_name), 1);
        if isempty(probability_column)
            continue;
        end
        hit = classification(:, probability_column) >= ...
            config.artifact_thresholds(class_index);
        remove_mask = remove_mask | hit;
        for ic_index = find(hit(:))'
            if strlength(removal_reason(ic_index)) == 0
                removal_reason(ic_index) = class_name;
            else
                removal_reason(ic_index) = removal_reason(ic_index) + ...
                    "+" + class_name;
            end
        end
    end

    report.RemovedICs = find(remove_mask)';
    EEG.etc.ic_classification = EEG_for_label.etc.ic_classification;
    save_iclabel_report(EEG_for_label, classification, classes, ...
        remove_mask, removal_reason, qc_folder, safe_base);
    if ~isempty(report.RemovedICs)
        EEG = pop_subcomp(EEG, report.RemovedICs, 0);
    end
end

function save_iclabel_report(EEG, classification, classes, ...
        remove_mask, removal_reason, qc_folder, safe_base)
    n_components = size(classification, 1);
    report_table = table((1:n_components)', ...
        'VariableNames', {'IC'});
    for class_index = 1:numel(classes)
        variable_name = matlab.lang.makeValidName(classes(class_index));
        report_table.(variable_name) = ...
            classification(:, class_index) * 100;
    end
    report_table.Removed = remove_mask;
    report_table.RemovalReason = removal_reason;
    writetable(report_table, ...
        fullfile(qc_folder, safe_base + "_ICLabel.csv"), ...
        'Encoding', 'UTF-8');

    try
        component_activity = eeg_getdatact(EEG, ...
            'component', 1:n_components);
        [spectrum, spectrum_freqs] = spectopo(component_activity, ...
            0, EEG.srate, 'plot', 'off');
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [50 50 1600 1000]);
        cleanup = onCleanup(@() close(fig));
        for ic_index = 1:n_components
            subplot(4, 4, 2 * ic_index - 1);
            topoplot(EEG.icawinv(:, ic_index), ...
                EEG.chanlocs(EEG.icachansind), ...
                'electrodes', 'on', 'style', 'both', 'headrad', 0.5);
            [maximum_probability, maximum_index] = ...
                max(classification(ic_index, :));
            rejection_mark = "";
            if remove_mask(ic_index)
                rejection_mark = " [除去]";
            end
            title(sprintf('IC %d: %s %.1f%%%s', ic_index, ...
                classes(maximum_index), maximum_probability * 100, ...
                rejection_mark), 'Interpreter', 'none');

            subplot(4, 4, 2 * ic_index);
            plot(spectrum_freqs, spectrum(ic_index, :), ...
                'LineWidth', 1.5);
            xlim([1 200]);
            grid on;
            title(sprintf('IC %d Spectrum', ic_index));
            if ic_index >= 7
                xlabel('Frequency (Hz)');
            end
        end
        exportgraphics(fig, ...
            fullfile(qc_folder, safe_base + "_ICLabel.png"), ...
            'Resolution', 180);
        clear cleanup;
    catch ME
        warning('run_set_cwt_band_timecourse:ICLabelFigureFailed', ...
            'ICLabel画像の保存に失敗しました: %s', ME.message);
    end
end

function [ersp, freqs, valid_event_count] = ...
        compute_low_frequency_ersp(EEG, channel_index, ...
            event_latencies, config)
    % timefreqは「時間点 x 試行」を期待する。行ベクトルを渡すと
    % 1時間点 x 多数試行と解釈され、巨大配列を確保してしまう。
    signal = double(reshape(EEG.data(channel_index, :, :), [], 1));
    requested_freqs = config.low_frequencies_hz;
    longest_wavelet_samples = ceil( ...
        config.wavelet_cycles * EEG.srate / min(requested_freqs));
    requested_time_points = floor( ...
        (EEG.xmax - EEG.xmin) * 1000 / config.output_step_ms) + 1;
    maximum_time_points = EEG.pnts - longest_wavelet_samples;
    n_timesout = min(requested_time_points, maximum_time_points);
    if n_timesout < 10
        error('run_set_cwt_band_timecourse:RecordingTooShort', ...
            '0.5 Hz・3サイクルCWTに対して連続記録が短すぎます。');
    end

    [tf, freqs, continuous_times_ms] = timefreq(signal, EEG.srate, ...
        'cycles', config.wavelet_cycles, ...
        'freqs', requested_freqs, ...
        'ntimesout', n_timesout, ...
        'tlimits', [EEG.xmin EEG.xmax] * 1000, ...
        'padratio', config.padratio, ...
        'verbose', 'off');
    tf = squeeze(tf);
    if size(tf, 1) ~= numel(freqs)
        tf = tf.';
    end
    power_values = abs(tf).^2;

    display_times = config.common_times_ms;
    baseline_times = config.baseline_window_ms(1):...
        config.output_step_ms:config.baseline_window_ms(2);
    n_events = numel(event_latencies);
    event_power = nan(numel(freqs), numel(display_times), n_events);
    baseline_power = nan(numel(freqs), numel(baseline_times), n_events);
    valid = false(1, n_events);
    event_times_ms = EEG.xmin * 1000 + ...
        (double(event_latencies(:)') - 1) / EEG.srate * 1000;

    for event_index = 1:n_events
        display_query = event_times_ms(event_index) + display_times;
        baseline_query = event_times_ms(event_index) + baseline_times;
        all_queries = [display_query, baseline_query];
        if min(all_queries) < continuous_times_ms(1) || ...
                max(all_queries) > continuous_times_ms(end)
            continue;
        end
        event_power(:, :, event_index) = interp1( ...
            continuous_times_ms(:), power_values.', ...
            display_query(:), 'linear').';
        baseline_power(:, :, event_index) = interp1( ...
            continuous_times_ms(:), power_values.', ...
            baseline_query(:), 'linear').';
        valid(event_index) = true;
    end

    valid_event_count = sum(valid);
    if valid_event_count < 2
        error('run_set_cwt_band_timecourse:TooFewLowFrequencyEvents', ...
            '低周波CWTの端部条件を満たすSaccadeが不足しています。');
    end
    event_power = event_power(:, :, valid);
    baseline_power = baseline_power(:, :, valid);
    mean_event_power = mean(event_power, 3, 'omitnan');
    baseline_reference = mean(reshape(baseline_power, ...
        numel(freqs), []), 2, 'omitnan');
    ersp = 10 * log10(max(mean_event_power, realmin) ./ ...
        max(baseline_reference, realmin));
    freqs = double(freqs(:));
end

function [ersp_common, freqs] = ...
        compute_high_frequency_ersp(EEG, channel_index, config)
    [ersp, ~, ~, times, freqs] = newtimef( ...
        EEG.data(channel_index, :, :), EEG.pnts, ...
        [EEG.xmin EEG.xmax] * 1000, EEG.srate, ...
        config.wavelet_cycles, ...
        'plotphase', 'off', 'plotersp', 'off', 'plotitc', 'off', ...
        'baseline', config.baseline_window_ms, ...
        'freqs', config.high_frequencies_hz, ...
        'timesout', config.high_timesout, ...
        'padratio', config.padratio, ...
        'verbose', 'off');
    ersp = double(ersp);
    freqs = double(freqs(:));
    ersp_common = interp1(double(times(:)), ersp.', ...
        config.common_times_ms(:), 'linear', NaN).';
end

function [band_timecourses, summary, timecourse_table] = ...
        calculate_band_timecourses(freqs, ersp, times, config)
    n_bands = numel(config.band_names);
    n_times = numel(times);
    band_timecourses = nan(n_bands, n_times);
    bin_count = zeros(n_bands, 1);
    excluded_line_bin_count = zeros(n_bands, 1);
    coverage = strings(n_bands, 1);
    actual_minimum = nan(n_bands, 1);
    actual_maximum = nan(n_bands, 1);
    mean_ersp = nan(n_bands, 1);

    for band_index = 1:n_bands
        low_hz = config.band_ranges(band_index, 1);
        high_hz = config.band_ranges(band_index, 2);
        original_mask = freqs >= low_hz & freqs <= high_hz;
        band_mask = original_mask;
        if config.exclude_line_noise_bins
            for line_frequency = config.line_frequencies_hz
                band_mask = band_mask & ...
                    abs(freqs - line_frequency) > ...
                    config.line_exclusion_half_width_hz;
            end
        end
        excluded_line_bin_count(band_index) = ...
            sum(original_mask) - sum(band_mask);
        bin_count(band_index) = sum(band_mask);
        if ~any(band_mask)
            coverage(band_index) = "なし";
            continue;
        end
        actual_minimum(band_index) = min(freqs(band_mask));
        actual_maximum(band_index) = max(freqs(band_mask));
        band_timecourses(band_index, :) = ...
            mean(ersp(band_mask, :), 1, 'omitnan');
        mean_ersp(band_index) = mean( ...
            band_timecourses(band_index, :), 'omitnan');
        if min(freqs) <= low_hz && max(freqs) >= high_hz
            coverage(band_index) = "全域";
        else
            coverage(band_index) = "一部";
        end
    end

    summary = table(config.band_names, config.band_labels, ...
        config.band_ranges(:, 1), config.band_ranges(:, 2), ...
        coverage, config.band_interpretation, bin_count, ...
        excluded_line_bin_count, actual_minimum, actual_maximum, ...
        mean_ersp, ...
        'VariableNames', {'Band', 'BandLabelJapanese', ...
            'DefinedMinimumHz', 'DefinedMaximumHz', ...
            'Coverage', 'Interpretation', 'FrequencyBinCount', ...
            'ExcludedLineNoiseBinCount', 'ActualMinimumHz', ...
            'ActualMaximumHz', 'MeanERSPdB'});

    time_ms = repmat(times(:), n_bands, 1);
    band = repelem(config.band_names, n_times);
    band_label = repelem(config.band_labels, n_times);
    interpretation = repelem(config.band_interpretation, n_times);
    mean_ersp_timecourse = reshape(band_timecourses.', [], 1);
    timecourse_table = table(time_ms, band, band_label, ...
        interpretation, mean_ersp_timecourse, ...
        'VariableNames', {'TimeMs', 'Band', 'BandLabelJapanese', ...
            'Interpretation', 'MeanERSPdB'});
end

function manifest = create_green_difference_outputs( ...
        folder_results, output_root, config)
    manifest = table();
    if isempty(folder_results)
        return;
    end

    n_results = numel(folder_results);
    is_noled = false(n_results, 1);
    for result_index = 1:n_results
        is_noled(result_index) = contains( ...
            string(folder_results{result_index}.BaseName), "緑なし");
    end
    noled_indices = find(is_noled);
    if isempty(noled_indices)
        fprintf('  緑なしSETがないため、ERSP差分出力はありません。\n');
        return;
    end

    condition_indices = find(~is_noled);
    fprintf('\n  --- 緑なしERSP差分を作成 ---\n');
    for condition_index = condition_indices(:)'
        condition = folder_results{condition_index};
        [noled_index, pairing_method] = select_noled_reference( ...
            folder_results, noled_indices, condition_index);
        if isempty(noled_index)
            message = "対応する緑なしSETを一意に決定できません。";
            warning('run_set_cwt_band_timecourse:NoUniqueNoLEDReference', ...
                '%s: %s', condition.BaseName, message);
            for channel_index = 1:numel(condition.ChannelNames)
                row = make_difference_manifest_row(condition, struct(), ...
                    condition.ChannelNames(channel_index), "失敗", ...
                    message, "", "", "", "", "");
                manifest = append_table_row(manifest, row);
            end
            continue;
        end

        noled = folder_results{noled_index};
        fprintf('  %s - %s (%s)\n', condition.BaseName, ...
            noled.BaseName, pairing_method);
        for channel_index = 1:numel(condition.ChannelNames)
            channel_name = condition.ChannelNames(channel_index);
            try
                [cwt_mat, cwt_png, band_mat, band_png] = ...
                    write_one_green_difference(condition, noled, ...
                        channel_index, output_root, config);
                row = make_difference_manifest_row(condition, noled, ...
                    channel_name, "完了", "", pairing_method, ...
                    cwt_mat, cwt_png, band_mat, band_png);
            catch ME
                warning('run_set_cwt_band_timecourse:DifferenceFailed', ...
                    '%s %sの緑なし差分に失敗しました: %s', ...
                    condition.BaseName, channel_name, ME.message);
                row = make_difference_manifest_row(condition, noled, ...
                    channel_name, "失敗", string(ME.message), ...
                    pairing_method, "", "", "", "");
            end
            manifest = append_table_row(manifest, row);
        end
    end
end

function [noled_index, method] = select_noled_reference( ...
        folder_results, noled_indices, condition_index)
    noled_index = [];
    method = "";
    if isscalar(noled_indices)
        noled_index = noled_indices;
        method = "単一緑なしSET";
        return;
    end

    condition_name = folder_results{condition_index}.BaseName;
    exact_condition_key = make_set_pairing_key(condition_name, false);
    exact_noled_keys = strings(numel(noled_indices), 1);
    for index = 1:numel(noled_indices)
        exact_noled_keys(index) = make_set_pairing_key( ...
            folder_results{noled_indices(index)}.BaseName, false);
    end
    exact_hits = noled_indices(exact_noled_keys == exact_condition_key);
    if isscalar(exact_hits)
        noled_index = exact_hits;
        method = "被験者・測定日一致";
        return;
    end

    relaxed_condition_key = make_set_pairing_key(condition_name, true);
    relaxed_noled_keys = strings(numel(noled_indices), 1);
    for index = 1:numel(noled_indices)
        relaxed_noled_keys(index) = make_set_pairing_key( ...
            folder_results{noled_indices(index)}.BaseName, true);
    end
    relaxed_hits = noled_indices( ...
        relaxed_noled_keys == relaxed_condition_key);
    if isscalar(relaxed_hits)
        noled_index = relaxed_hits;
        method = "被験者一致（測定日差を許容）";
    end
end

function key = make_set_pairing_key(base_name, ignore_measurement_date)
    key = regexprep(string(base_name), '緑なし', '');
    key = regexprep(key, '(^|_)[0-9]+[Oo](_|$)', '$1$2');
    if ignore_measurement_date
        key = regexprep(key, '\s*[\(\uff08][0-9]+[\)\uff09]', '');
        key = regexprep(key, '(^|_)[0-9]{4}(_|$)', '$1$2');
    end
    key = regexprep(key, '[_\s-]+', '_');
    key = regexprep(key, '^_+|_+$', '');
    key = lower(key);
end

function [output_cwt_mat, output_cwt_png, ...
        output_band_mat, output_band_png] = ...
        write_one_green_difference(condition, noled, channel_index, ...
            output_root, config)
    condition_data = load(condition.CWTFiles(channel_index), ...
        'ersp', 'freqs', 'times');
    noled_data = load(noled.CWTFiles(channel_index), ...
        'ersp', 'freqs', 'times');
    validate_difference_inputs(condition_data, noled_data, ...
        condition.CWTFiles(channel_index), noled.CWTFiles(channel_index));

    ersp = double(condition_data.ersp) - double(noled_data.ersp);
    freqs = double(condition_data.freqs(:));
    times = double(condition_data.times(:)');
    [band_timecourses, band_summary, timecourse_table] = ...
        calculate_band_timecourses(freqs, ersp, times, config);

    condition_source_set = string(condition.SourceSET);
    green_noled_source_set = string(noled.SourceSET);
    analysis_config = config;
    difference_metadata = struct( ...
        'quantity', "ERSP difference in dB", ...
        'formula', "condition ERSP - green-noLED ERSP", ...
        'condition_source_set', condition_source_set, ...
        'green_noled_source_set', green_noled_source_set, ...
        'note', "両条件をそれぞれベースライン補正した" + ...
                "ERSP(dB)間の差です。");

    channel_name = condition.ChannelNames(channel_index);
    channel_output = fullfile(output_root, channel_name);
    if ~isfolder(channel_output)
        mkdir(channel_output);
    end
    safe_base = sanitize_filename(condition.BaseName);
    cwt_stem = safe_base + "_" + channel_name + ...
        "_CWT再計算_緑なし差分";
    band_stem = safe_base + "_" + channel_name + ...
        "_SET再計算_緑なし差分_band_timecourse";
    output_cwt_mat = fullfile(channel_output, cwt_stem + ".mat");
    output_cwt_png = fullfile(channel_output, cwt_stem + ".png");
    output_band_mat = fullfile(channel_output, band_stem + ".mat");
    output_band_png = fullfile(channel_output, band_stem + ".png");
    output_timecourse_csv = fullfile(channel_output, band_stem + ".csv");
    output_summary_csv = fullfile(channel_output, ...
        band_stem + "_summary.csv");

    low_mask = freqs <= max(config.low_frequencies_hz);
    low_ersp = ersp(low_mask, :);
    low_freqs = freqs(low_mask);
    high_ersp = ersp(~low_mask, :);
    high_freqs = freqs(~low_mask);
    timecourse_table.Properties.VariableNames{'MeanERSPdB'} = ...
        'MeanERSPDifferenceDb';
    band_summary.Properties.VariableNames{'MeanERSPdB'} = ...
        'MeanERSPDifferenceDb';
    save(output_cwt_mat, 'condition_source_set', ...
        'green_noled_source_set', 'difference_metadata', ...
        'analysis_config', 'ersp', 'freqs', 'times', 'low_ersp', ...
        'low_freqs', 'high_ersp', 'high_freqs', '-v7');
    save(output_band_mat, 'condition_source_set', ...
        'green_noled_source_set', 'difference_metadata', ...
        'analysis_config', 'times', 'band_timecourses', ...
        'band_summary', '-v7');

    writetable(timecourse_table, output_timecourse_csv, ...
        'Encoding', 'UTF-8');
    writetable(band_summary, output_summary_csv, ...
        'Encoding', 'UTF-8');

    source_label = condition.BaseName + " - " + noled.BaseName;
    create_band_figure(times, band_timecourses, band_summary, ...
        source_label, channel_name, output_band_png, config, true);
    create_cwt_figure(times, freqs, ersp, source_label, ...
        channel_name, output_cwt_png, config, true);
end

function validate_difference_inputs(condition_data, noled_data, ...
        condition_file, noled_file)
    required_fields = {'ersp', 'freqs', 'times'};
    for field_index = 1:numel(required_fields)
        field_name = required_fields{field_index};
        if ~isfield(condition_data, field_name) || ...
                ~isfield(noled_data, field_name)
            error('run_set_cwt_band_timecourse:MissingDifferenceField', ...
                '%sまたは%sに%sがありません。', ...
                condition_file, noled_file, field_name);
        end
    end
    if ~isequal(size(condition_data.ersp), size(noled_data.ersp))
        error('run_set_cwt_band_timecourse:DifferenceSizeMismatch', ...
            '条件と緑なしのERSPサイズが一致しません。');
    end
    condition_freqs = double(condition_data.freqs(:));
    noled_freqs = double(noled_data.freqs(:));
    condition_times = double(condition_data.times(:));
    noled_times = double(noled_data.times(:));
    if numel(condition_freqs) ~= numel(noled_freqs) || ...
            any(abs(condition_freqs - noled_freqs) > 1e-9)
        error('run_set_cwt_band_timecourse:DifferenceFrequencyMismatch', ...
            '条件と緑なしの周波数軸が一致しません。');
    end
    if numel(condition_times) ~= numel(noled_times) || ...
            any(abs(condition_times - noled_times) > 1e-9)
        error('run_set_cwt_band_timecourse:DifferenceTimeMismatch', ...
            '条件と緑なしの時間軸が一致しません。');
    end
end

function row = make_difference_manifest_row(condition, noled, ...
        channel_name, status, message, pairing_method, cwt_mat, ...
        cwt_png, band_mat, band_png)
    if isempty(fieldnames(noled))
        noled_source = "";
    else
        noled_source = string(noled.SourceSET);
    end
    row = table(string(condition.SourceSET), noled_source, ...
        string(channel_name), string(status), string(message), ...
        string(pairing_method), string(cwt_mat), string(cwt_png), ...
        string(band_mat), string(band_png), ...
        'VariableNames', {'ConditionSET', 'GreenNoLEDSET', 'Channel', ...
            'Status', 'Message', 'PairingMethod', 'OutputCWTMAT', ...
            'OutputCWTPNG', 'OutputBandMAT', 'OutputBandPNG'});
end

function output = append_table_row(output, row)
    if isempty(output)
        output = row;
    else
        output = [output; row];
    end
end

function create_band_figure(times, band_timecourses, summary, ...
        source_base, channel_name, output_png, config, is_difference)
    if nargin < 8
        is_difference = false;
    end
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [80 40 1350 1200]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 4, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    finite_values = band_timecourses(isfinite(band_timecourses));
    if isempty(finite_values)
        common_ylim = [-1 1];
    else
        common_ylim = [min(finite_values), max(finite_values)];
        if common_ylim(1) == common_ylim(2)
            common_ylim = common_ylim + [-1 1];
        else
            margin = max(0.1, 0.08 * diff(common_ylim));
            common_ylim = common_ylim + [-margin margin];
        end
    end

    for band_index = 1:height(summary)
        ax = nexttile(layout);
        if summary.Coverage(band_index) ~= "なし"
            line_color = [0.05 0.35 0.70];
            if is_difference
                line_color = [0.75 0.20 0.10];
            end
            plot(ax, times, band_timecourses(band_index, :), ...
                'LineWidth', 1.7, 'Color', line_color);
            hold(ax, 'on');
            yline(ax, 0, '--k', 'LineWidth', 0.8);
            ylim(ax, common_ylim);
        else
            ylim(ax, [-1 1]);
            text(ax, mean(config.display_window_ms), 0, ...
                'データなし', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', [0.75 0.15 0.10], 'FontWeight', 'bold');
            hold(ax, 'on');
        end
        xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1);
        grid(ax, 'on');
        xlim(ax, config.display_window_ms);
        xlabel(ax, '時間 (ms)');
        if is_difference
            ylabel(ax, '平均ERSP差分 (dB)');
        else
            ylabel(ax, '平均ERSP (dB)');
        end
        title(ax, sprintf('%s (%.1f–%.1f Hz) [%s]', ...
            summary.BandLabelJapanese(band_index), ...
            summary.DefinedMinimumHz(band_index), ...
            summary.DefinedMaximumHz(band_index), ...
            summary.Coverage(band_index)), ...
            'Interpreter', 'none');
        set(ax, 'FontSize', 10, 'Box', 'on');
    end

    if is_difference
        quantity_label = '帯域平均ERSP差分（条件 - 緑なし）';
    else
        quantity_label = '帯域平均ERSP時系列';
    end
    title(layout, sprintf( ...
        ['%s_%s / %s: %s ' ...
         '(%.0f～%.0f ms, 0 ms = サッケード開始)'], ...
        source_base, channel_name, channel_name, quantity_label, ...
        config.display_window_ms(1), config.display_window_ms(2)), ...
        'Interpreter', 'none', 'FontWeight', 'bold');
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function create_cwt_figure(times, freqs, ersp, source_base, ...
        channel_name, output_png, config, is_difference)
    if nargin < 8
        is_difference = false;
    end
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [80 80 1250 700]);
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);
    surface(ax, times, freqs, zeros(size(ersp)), ersp, ...
        'EdgeColor', 'none', 'FaceColor', 'flat');
    view(ax, 2);
    axis(ax, 'tight');
    xlim(ax, config.display_window_ms);
    ylim(ax, [min(freqs) max(freqs)]);
    colormap(ax, turbo);
    colorbar_handle = colorbar(ax);
    if is_difference
        colorbar_handle.Label.String = 'ERSP差分 (dB)';
    else
        colorbar_handle.Label.String = 'ERSP (dB)';
    end
    finite_values = ersp(isfinite(ersp));
    if ~isempty(finite_values)
        symmetric_limit = prctile(abs(finite_values), 98);
        if isfinite(symmetric_limit) && symmetric_limit > 0
            clim(ax, [-symmetric_limit symmetric_limit]);
        end
    end
    hold(ax, 'on');
    xline(ax, 0, '--m', 'LineWidth', 1.3);
    xlabel(ax, '時間 (ms)');
    ylabel(ax, '周波数 (Hz)');
    if is_difference
        quantity_label = 'ERSP差分（条件 - 緑なし）';
    else
        quantity_label = 'CWT再計算ERSP';
    end
    title(ax, sprintf( ...
        ['%s_%s / %s: %s ' ...
         '(%.0f～%.0f ms, 0 ms = サッケード開始)'], ...
        source_base, channel_name, channel_name, quantity_label, ...
        config.display_window_ms(1), config.display_window_ms(2)), ...
        'Interpreter', 'none', 'FontWeight', 'bold');
    set(ax, 'FontSize', 11, 'Box', 'on');
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function latencies = get_event_latencies(EEG, event_type)
    types = event_types_as_strings(EEG.event);
    mask = strcmpi(types, event_type);
    latencies = double([EEG.event(mask).latency]);
    latencies = sort(latencies);
end

function types = event_types_as_strings(events)
    types = strings(1, numel(events));
    for index = 1:numel(events)
        types(index) = string(events(index).type);
    end
end

function row = make_manifest_row(manual_folder, set_path, ...
        output_root, status, message, result)
    defaults = struct( ...
        'RawSaccadeEvents', NaN, ...
        'UniqueSaccadeEvents', NaN, ...
        'DuplicatesRemoved', NaN, ...
        'HighEpochCount', NaN, ...
        'LowValidEventCount', NaN, ...
        'RemovedICs', "");
    names = fieldnames(defaults);
    for index = 1:numel(names)
        if ~isfield(result, names{index})
            result.(names{index}) = defaults.(names{index});
        end
    end
    row = table(string(manual_folder), string(set_path), ...
        string(output_root), string(status), string(message), ...
        result.RawSaccadeEvents, result.UniqueSaccadeEvents, ...
        result.DuplicatesRemoved, result.HighEpochCount, ...
        result.LowValidEventCount, string(result.RemovedICs), ...
        'VariableNames', {'ManualSaccadeFolder', 'SourceSET', ...
            'OutputRoot', 'Status', 'Message', 'RawSaccadeEvents', ...
            'UniqueSaccadeEvents', 'DuplicatesRemoved', ...
            'HighEpochCount', 'LowValidEventCount', 'RemovedICs'});
end

function safe_name = sanitize_filename(name)
    safe_name = regexprep(string(name), '[<>:"/|?*]', '_');
    safe_name = strrep(safe_name, '\', '_');
end
