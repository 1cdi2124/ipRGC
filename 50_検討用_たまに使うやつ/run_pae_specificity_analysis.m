function master_manifest = run_pae_specificity_analysis( ...
        input_root, analysis_channels)
%RUN_PAE_SPECIFICITY_ANALYSIS PAE候補の特異性を検討する探索解析。
%
% 実行例:
%   run_pae_specificity_analysis
%   run_pae_specificity_analysis("C:\...\03_手動サッケード")
%   run_pae_specificity_analysis("C:\...\03_手動サッケード", "Oz")
%
% 各03_手動サッケード内のSETを同一被験者の条件として扱う。
% 0Oは緑LED常灯、緑なしは無灯の眼球運動対照であり、区別する。
% 全条件へ同一フィルタ・平均基準を適用し、連結データでICAを1回だけ学習する。
% 本感度分析では、ICLabel Brain確率50%以上かつ
% |IC-EOG相関|<0.40のICだけを残し、同じIC選択を全条件へ適用する。
% 条件別のTotal/Induced/Evoked ERSP、ITC、試行EOG指標を保存し、
% Alpha以上についてEOG共変量調整後の条件差も出力する。
% Delta/Thetaも可視化するが、反復サッケード重畳の探索的指標である。

    config = default_config();
    if nargin >= 2 && ~isempty(analysis_channels)
        analysis_channels = string(analysis_channels(:)');
        if ~all(ismember(lower(analysis_channels), ...
                lower(config.target_channels)))
            error('run_pae_specificity_analysis:InvalidOutputChannel', ...
                '解析チャンネルは次から選んでください: %s', ...
                strjoin(config.target_channels, ', '));
        end
        config.analysis_channels = analysis_channels;
    end
    if nargin < 1 || strlength(string(input_root)) == 0
        selected = uigetdir('', ...
            ['実験ルート、被験者フォルダ、または' ...
             '03_手動サッケードを選択してください']);
        if isequal(selected, 0)
            error('run_pae_specificity_analysis:Cancelled', ...
                '処理を中断しました。');
        end
        input_root = string(selected);
    else
        input_root = string(input_root);
    end
    if ~isfolder(input_root)
        error('run_pae_specificity_analysis:MissingInputFolder', ...
            '入力フォルダが見つかりません: %s', input_root);
    end

    initialize_eeglab(config);
    manual_folders = find_manual_saccade_folders(input_root);
    if isempty(manual_folders)
        error('run_pae_specificity_analysis:NoInputFolder', ...
            '03_手動サッケードが見つかりません: %s', input_root);
    end

    master_manifest = table();
    fprintf('\n=== PAE特異性解析（Brain 50%% 共通ICA）を開始 ===\n');
    fprintf('検索ルート: %s\n', input_root);
    fprintf('対象フォルダ数: %d\n', numel(manual_folders));
    fprintf('解析フィルタ: %.1f–%.0f Hz / 平均基準（EOG除外）\n', ...
        config.analysis_highpass_hz, config.analysis_lowpass_hz);

    for folder_index = 1:numel(manual_folders)
        manual_folder = manual_folders(folder_index);
        subject_folder = string(fileparts(manual_folder));
        output_root = fullfile(subject_folder, ...
            config.output_relative_path);
        ensure_folder(output_root);
        fprintf('\n--- [%d/%d] %s ---\n', folder_index, ...
            numel(manual_folders), manual_folder);

        try
            datasets = load_and_preprocess_conditions( ...
                manual_folder, config);
            [datasets, joint_ica_report] = run_joint_ica_cleaning( ...
                datasets, output_root, config);
            [folder_manifest, analysis_results] = ...
                analyze_all_conditions(datasets, output_root, ...
                    joint_ica_report, config);
            create_all_channel_overlays(analysis_results, ...
                output_root, config);
            create_pae_specificity_outputs(analysis_results, ...
                output_root, config);
            folder_status = "完了";
            folder_message = "";
        catch ME
            warning('run_pae_specificity_analysis:FolderFailed', ...
                '%s の処理に失敗しました: %s', ...
                manual_folder, ME.message);
            folder_manifest = table();
            folder_status = "失敗";
            folder_message = string(getReport(ME, 'extended', ...
                'hyperlinks', 'off'));
        end

        folder_row = table(manual_folder, output_root, ...
            folder_status, folder_message, ...
            'VariableNames', {'ManualSaccadeFolder', 'OutputRoot', ...
                'Status', 'Message'});
        if isempty(master_manifest)
            master_manifest = folder_row;
        else
            master_manifest = [master_manifest; folder_row]; %#ok<AGROW>
        end
        if ~isempty(folder_manifest)
            writetable(folder_manifest, ...
                fullfile(output_root, 'PAE特異性解析_manifest.csv'), ...
                'Encoding', 'UTF-8');
        end
    end

    fprintf('\n=== PAE特異性解析（Brain 50%% 共通ICA）が終了しました ===\n');
    fprintf('完了フォルダ: %d / 失敗フォルダ: %d\n', ...
        sum(master_manifest.Status == "完了"), ...
        sum(master_manifest.Status == "失敗"));
    fprintf(['注意: Delta/Thetaは削除せず可視化していますが、' ...
        '反復サッケード重畳を含む探索的指標です。\n']);
end

function config = default_config()
    config.output_relative_path = "07_PAE特異性解析_Brain50";
    config.target_channels = ["F3", "F4", "Fz", "O1", ...
        "O2", "Oz", "PO7", "PO8"];
    config.analysis_channels = config.target_channels;
    config.eog_channel_patterns = ["EOG", "HEOG", "VEOG"];
    config.event_type = "Saccade";
    config.event_duplicate_tolerance_samples = 1;

    config.analysis_highpass_hz = 0.2;
    config.analysis_lowpass_hz = 220;
    config.apply_average_reference = true;
    config.ica_highpass_hz = 1;
    config.ica_lowpass_hz = 100;
    config.ica_training_srate = 250;
    config.ica_rank = 7;
    config.ica_seed = 3;
    config.brain_keep_threshold = 0.50;
    config.eog_correlation_threshold = 0.40;

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
        "探索的（反復サッケード重畳）"; ...
        "探索的（反復サッケード重畳）"; ...
        "可視化"; "可視化"; "可視化"; "可視化"; ...
        "探索的"; "探索的"];
    config.confidence_multiplier = 1.96;
    config.specificity_band_indices = [3 4 5 6 7 8];
    config.specificity_windows_ms = [ ...
          0, 150; ... % Alpha
          0, 150; ... % Beta
          0,  50; ... % Low gamma
          0, 150; ... % High gamma 1
          0, 150; ... % High gamma 2
          0, 150];    % High gamma 3
    config.eog_metric_window_ms = [-100 200];

    config.eeglab_fallback = ...
        "C:\Users\tatsuya\Downloads\Apps\eeglab2025.1.0";
end

function initialize_eeglab(config)
    if exist('eeglab', 'file') ~= 2 && isfolder(config.eeglab_fallback)
        addpath(config.eeglab_fallback);
    end
    if exist('eeglab', 'file') ~= 2
        error('run_joint_ica_cwt_band_analysis:EEGLABNotFound', ...
            'EEGLABがMATLABパスにありません。');
    end
    if exist('timefreq', 'file') ~= 2 || ...
            exist('pop_eegfiltnew', 'file') ~= 2
        evalc('eeglab nogui');
    elseif exist('pop_iclabel', 'file') ~= 2
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

function datasets = load_and_preprocess_conditions(manual_folder, config)
    set_list = dir(fullfile(manual_folder, '*.set'));
    if numel(set_list) < 2
        error('run_joint_ica_cwt_band_analysis:TooFewConditions', ...
            '共通ICAには2条件以上のSETが必要です: %s', manual_folder);
    end
    datasets = repmat(struct(), numel(set_list), 1);
    reference_labels = strings(0, 1);
    fprintf('SET数: %d\n', numel(set_list));

    for file_index = 1:numel(set_list)
        set_path = string(fullfile(set_list(file_index).folder, ...
            set_list(file_index).name));
        [~, base_name] = fileparts(set_path);
        fprintf('  読込・共通前処理: %s\n', base_name);
        EEG = load_set_with_fdt_repair(set_path);
        EEG = clear_ica_fields(EEG);
        if EEG.trials ~= 1
            error('run_joint_ica_cwt_band_analysis:NotContinuous', ...
                '%s は連続データではありません。', set_path);
        end
        [EEG, raw_count, removed_count] = ...
            deduplicate_saccade_events(EEG, config.event_type, ...
                config.event_duplicate_tolerance_samples);
        EEG = remove_nonfinite_segments(EEG);

        labels = string({EEG.chanlocs.labels});
        if isempty(reference_labels)
            reference_labels = labels;
        elseif ~isequal(lower(labels), lower(reference_labels))
            error('run_joint_ica_cwt_band_analysis:ChannelMismatch', ...
                '条件間でチャンネル順が一致しません: %s', set_path);
        end
        target_indices = find_target_channels( ...
            EEG, config.target_channels);
        eog_index = find_eog_channel(EEG, config.eog_channel_patterns);

        EEG = pop_eegfiltnew(EEG, ...
            'locutoff', config.analysis_highpass_hz, ...
            'hicutoff', config.analysis_lowpass_hz, ...
            'channels', 1:EEG.nbchan, 'plotfreqz', 0);
        if config.apply_average_reference
            EEG = pop_reref(EEG, [], 'exclude', eog_index);
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');

        [condition_label, condition_value, is_green] = ...
            parse_condition_name(base_name);
        datasets(file_index).EEG = EEG;
        datasets(file_index).SourceSET = set_path;
        datasets(file_index).BaseName = string(base_name);
        datasets(file_index).ConditionLabel = condition_label;
        datasets(file_index).ConditionValue = condition_value;
        datasets(file_index).IsGreen = is_green;
        datasets(file_index).TargetIndices = target_indices;
        datasets(file_index).EOGIndex = eog_index;
        datasets(file_index).RawSaccadeEvents = raw_count;
        datasets(file_index).DuplicatesRemoved = removed_count;
    end

    sort_values = [datasets.ConditionValue]';
    sort_values([datasets.IsGreen]') = Inf;
    [~, order] = sort(sort_values);
    datasets = datasets(order);
end

function [datasets, report] = run_joint_ica_cleaning( ...
        datasets, output_root, config)
    qc_folder = fullfile(output_root, 'ICA_QC');
    ensure_folder(qc_folder);
    [merged, segment_starts, segment_stops] = ...
        concatenate_condition_data(datasets);
    target_indices = datasets(1).TargetIndices;
    eog_index = datasets(1).EOGIndex;

    try
        if ~isfield(merged.chanlocs, 'X') || ...
                isempty(merged.chanlocs(target_indices(1)).X)
            lookup_file = resolve_standard_lookup_file();
            merged = pop_chanedit(merged, ...
                'lookup', char(lookup_file));
        end
    catch ME
        error('run_joint_ica_cwt_band_analysis:ChannelLookup', ...
            ['ICLabelに必要な標準電極座標を補完できません。' ...
            'EEGLAB/DIPFITのsupport fileを確認してください: %s'], ...
            ME.message);
    end

    training = merged;
    training = pop_eegfiltnew(training, ...
        'locutoff', config.ica_highpass_hz, ...
        'hicutoff', config.ica_lowpass_hz, ...
        'channels', 1:training.nbchan, 'plotfreqz', 0);
    if abs(training.srate - config.ica_training_srate) > 1e-6
        training = pop_resample(training, config.ica_training_srate);
    end

    rng(config.ica_seed, 'twister');
    data_rank = rank(double(training.data(target_indices, :)));
    rank_value = min([config.ica_rank, data_rank, ...
        numel(target_indices)]);
    ica_options = {'icatype', 'runica', 'chanind', target_indices, ...
        'extended', 1, 'rndreset', 'no'};
    ica_options = [ica_options, {'pca', rank_value}];
    fprintf('  共通ICA: %d EEG ch, rank=%d, 学習%.0f Hz\n', ...
        numel(target_indices), rank_value, training.srate);
    training = pop_runica(training, ica_options{:});

    merged.icaweights = training.icaweights;
    merged.icasphere = training.icasphere;
    merged.icawinv = training.icawinv;
    merged.icachansind = target_indices;
    merged.icaact = [];
    merged = eeg_checkset(merged, 'ica');

    if exist('pop_iclabel', 'file') == 2
        labeled = pop_iclabel(merged, 'default');
        classification = ...
            labeled.etc.ic_classification.ICLabel.classifications;
        classes = string( ...
            labeled.etc.ic_classification.ICLabel.classes);
        merged.etc.ic_classification = ...
            labeled.etc.ic_classification;
    else
        error('run_joint_ica_cwt_band_analysis:ICLabelMissing', ...
            ['明確なアーチファクトと要確認ICを分類するため、' ...
             'ICLabelが必要です。EEGLABへICLabelプラグインを' ...
             '追加してください。']);
    end

    training_ica = double(training.icaweights * ...
        training.icasphere * training.data(target_indices, :));
    training_eog = double(training.data(eog_index, :));
    eog_correlations = calculate_eog_correlations( ...
        training_ica, training_eog);
    [remove_mask, review_mask, removal_reason, review_reason] = ...
        choose_artifact_components(classification, classes, ...
        eog_correlations, config);
    removed_components = find(remove_mask(:))';
    review_components = find(review_mask(:))';
    kept_components = find(~remove_mask(:))';
    if isempty(removed_components)
        cleaned = merged;
    else
        cleaned = pop_subcomp(merged, removed_components, 0);
    end
    fprintf('  共通ICA採用IC: %s（要確認: %s）/ 除去IC: %s\n', ...
        mat2str(kept_components), mat2str(review_components), ...
        mat2str(removed_components));

    report = struct();
    report.Rank = rank_value;
    report.DataRank = data_rank;
    report.RequestedRank = config.ica_rank;
    report.Classes = classes;
    report.Classification = classification;
    report.EOGCorrelations = eog_correlations;
    report.KeptComponents = kept_components;
    report.ReviewComponents = review_components;
    report.RemovedComponents = removed_components;
    report.ReviewReason = review_reason;
    report.RemovalReason = removal_reason;
    report.AnalysisFilterHz = [config.analysis_highpass_hz, ...
        config.analysis_lowpass_hz];
    report.ICATrainingFilterHz = [config.ica_highpass_hz, ...
        config.ica_lowpass_hz];
    report.ICATrainingSrate = training.srate;
    save_joint_ica_qc(merged, training, report, qc_folder, config);

    for data_index = 1:numel(datasets)
        EEG = datasets(data_index).EEG;
        EEG.data = double(cleaned.data(:, ...
            segment_starts(data_index):segment_stops(data_index)));
        EEG.pnts = size(EEG.data, 2);
        EEG.xmax = EEG.xmin + (EEG.pnts - 1) / EEG.srate;
        EEG.times = linspace(EEG.xmin * 1000, ...
            EEG.xmax * 1000, EEG.pnts);
        EEG = clear_ica_fields(EEG);
        EEG = eeg_checkset(EEG, 'eventconsistency');
        datasets(data_index).EEG = EEG;
    end
end

function ensure_folder(folder_path)
    if ~isfolder(folder_path)
        [created, message] = mkdir(folder_path);
        if ~created
            error('run_joint_ica_cwt_band_analysis:CreateFolderFailed', ...
                '出力フォルダを作成できません: %s (%s)', ...
                folder_path, message);
        end
    end
end

function lookup_file = resolve_standard_lookup_file()
    eeglab_file = string(which('eeglab'));
    if strlength(eeglab_file) == 0
        error('run_joint_ica_cwt_band_analysis:EEGLABPathMissing', ...
            'eeglab.mの場所を取得できません。');
    end
    eeglab_root = string(fileparts(eeglab_file));
    candidates = [ ...
        fullfile(eeglab_root, 'plugins', 'dipfit', ...
            'standard_BESA', 'standard-10-5-cap385.elp'); ...
        fullfile(eeglab_root, 'functions', 'supportfiles', ...
            'Standard-10-5-Cap385.sfp'); ...
        fullfile(eeglab_root, 'functions', 'supportfiles', ...
            'Standard-10-5-Cap385_witheog.elp')];
    found = find(arrayfun(@isfile, candidates), 1);
    if isempty(found)
        error('run_joint_ica_cwt_band_analysis:LookupFileMissing', ...
            'EEGLAB内に標準10-5電極座標ファイルがありません。');
    end
    lookup_file = candidates(found);
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
                error('run_joint_ica_cwt_band_analysis:FDTNotFound', ...
                    '%s に対応するFDTを一意に決められません。', ...
                    set_path);
            end
        end
        EEG.data = eeg_getdatact(EEG);
    end
    EEG.data = double(EEG.data);
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

function EEG = clear_ica_fields(EEG)
    fields = {'icaact', 'icaweights', 'icasphere', ...
        'icawinv', 'icachansind'};
    for field_index = 1:numel(fields)
        EEG.(fields{field_index}) = [];
    end
    if isfield(EEG, 'etc') && ...
            isfield(EEG.etc, 'ic_classification')
        EEG.etc = rmfield(EEG.etc, 'ic_classification');
    end
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
        fprintf('    重複Saccadeイベントを%d個除去しました。\n', ...
            removed_count);
    end
end

function types = event_types_as_strings(events)
    types = strings(1, numel(events));
    for event_index = 1:numel(events)
        value = events(event_index).type;
        if iscell(value) && isscalar(value)
            value = value{1};
        end
        types(event_index) = string(value);
    end
end

function EEG = remove_nonfinite_segments(EEG)
    invalid = any(~isfinite(double(EEG.data)), 1);
    invalid = reshape(invalid, 1, []);
    if ~any(invalid)
        return;
    end
    valid = ~invalid;
    if sum(valid) < 2
        error('run_joint_ica_cwt_band_analysis:NonfiniteData', ...
            '有限値のデータ点が不足しています。');
    end
    x = 1:EEG.pnts;
    for channel_index = 1:EEG.nbchan
        signal = double(reshape(EEG.data(channel_index, :), 1, []));
        channel_valid = isfinite(signal) & valid;
        if sum(channel_valid) < 2
            error('run_joint_ica_cwt_band_analysis:NonfiniteChannel', ...
                'チャンネル%sの有限値が不足しています。', ...
                EEG.chanlocs(channel_index).labels);
        end
        signal(~channel_valid) = interp1(x(channel_valid), ...
            signal(channel_valid), x(~channel_valid), 'linear', ...
            'extrap');
        EEG.data(channel_index, :) = signal;
    end
    warning('run_joint_ica_cwt_band_analysis:InterpolatedSamples', ...
        '非有限値を含む%dサンプルを時間方向に補間しました。', ...
        sum(invalid));
end

function target_indices = find_target_channels(EEG, target_channels)
    labels = string({EEG.chanlocs.labels});
    target_indices = zeros(1, numel(target_channels));
    for channel_index = 1:numel(target_channels)
        found = find(strcmpi(labels, target_channels(channel_index)), 1);
        if isempty(found)
            error('run_joint_ica_cwt_band_analysis:MissingChannel', ...
                'チャンネル%sがSETにありません。', ...
                target_channels(channel_index));
        end
        target_indices(channel_index) = found;
    end
end

function eog_index = find_eog_channel(EEG, patterns)
    labels = string({EEG.chanlocs.labels});
    for pattern_index = 1:numel(patterns)
        found = find(strcmpi(labels, patterns(pattern_index)), 1);
        if ~isempty(found)
            eog_index = found;
            return;
        end
    end
    for pattern_index = 1:numel(patterns)
        found = find(contains(lower(labels), ...
            lower(patterns(pattern_index))), 1);
        if ~isempty(found)
            eog_index = found;
            return;
        end
    end
    error('run_joint_ica_cwt_band_analysis:EOGNotFound', ...
        'EOGチャンネルが見つかりません。候補: %s', ...
        strjoin(patterns, ', '));
end

function [label, value, is_green] = parse_condition_name(base_name)
    base_name = string(base_name);
    is_green = contains(base_name, "緑なし", 'IgnoreCase', true);
    if is_green
        label = "緑なし";
        value = NaN;
        return;
    end

    token = regexp(char(base_name), ...
        '(?:^|_)([0-9]+)[Oo](?:_|$)', 'tokens', 'once');
    if isempty(token)
        token = regexp(char(base_name), ...
            '(?:^|_)([0-9]+)Hz(?:_|$)', 'tokens', 'once', ...
            'ignorecase');
    end
    if isempty(token)
        label = base_name;
        value = NaN;
    else
        value = str2double(token{1});
        label = string(sprintf('%g Hz', value));
    end
end

function [merged, segment_starts, segment_stops] = ...
        concatenate_condition_data(datasets)
    n_datasets = numel(datasets);
    sample_counts = zeros(n_datasets, 1);
    reference = datasets(1).EEG;
    for data_index = 1:n_datasets
        EEG = datasets(data_index).EEG;
        if EEG.nbchan ~= reference.nbchan || ...
                abs(EEG.srate - reference.srate) > 1e-9
            error('run_joint_ica_cwt_band_analysis:DatasetMismatch', ...
                '条件間でチャンネル数またはサンプリング周波数が異なります。');
        end
        sample_counts(data_index) = EEG.pnts;
    end

    segment_starts = [1; cumsum(sample_counts(1:end-1)) + 1];
    segment_stops = cumsum(sample_counts);
    merged = reference;
    merged.data = zeros(reference.nbchan, sum(sample_counts), 'double');
    for data_index = 1:n_datasets
        merged.data(:, segment_starts(data_index):...
            segment_stops(data_index)) = double(datasets(data_index).EEG.data);
    end
    merged.trials = 1;
    merged.pnts = size(merged.data, 2);
    merged.xmin = 0;
    merged.xmax = (merged.pnts - 1) / merged.srate;
    merged.times = (0:merged.pnts - 1) / merged.srate * 1000;
    merged.event = struct([]);
    merged.urevent = struct([]);
    merged.epoch = [];
    merged.setname = '全条件連結_共通ICA';
    merged = clear_ica_fields(merged);
    merged = eeg_checkset(merged);
end

function correlations = calculate_eog_correlations(ica_activity, eog)
    n_components = size(ica_activity, 1);
    correlations = nan(n_components, 1);
    eog = double(eog(:));
    for component_index = 1:n_components
        activity = double(ica_activity(component_index, :))';
        valid = isfinite(activity) & isfinite(eog);
        if sum(valid) < 3 || std(activity(valid)) == 0 || ...
                std(eog(valid)) == 0
            continue;
        end
        coefficient = corrcoef(activity(valid), eog(valid));
        correlations(component_index) = coefficient(1, 2);
    end
end

function [remove_mask, review_mask, removal_reasons, ...
        review_reasons] = choose_artifact_components( ...
        classification, classes, eog_correlations, config)
    n_components = numel(eog_correlations);
    review_mask = false(n_components, 1);
    removal_reasons = strings(n_components, 1);
    review_reasons = strings(n_components, 1);
    if isempty(classification) || ...
            size(classification, 1) ~= n_components
        error('run_joint_ica_cwt_band_analysis:InvalidICLabelOutput', ...
            'ICLabelの成分数が共通ICAの成分数と一致しません。');
    end
    brain_column = find(strcmpi(classes, "Brain"), 1);
    if isempty(brain_column)
        error('run_joint_ica_cwt_band_analysis:BrainClassMissing', ...
            'ICLabel出力にBrainクラスがありません。');
    end

    brain_probability = classification(:, brain_column);
    brain_ok = isfinite(brain_probability) & ...
        brain_probability >= config.brain_keep_threshold;
    eog_ok = isfinite(eog_correlations) & ...
        abs(eog_correlations) < config.eog_correlation_threshold;
    keep_mask = brain_ok & eog_ok;
    remove_mask = ~keep_mask;

    for component_index = find(remove_mask(:))'
        if ~isfinite(brain_probability(component_index))
            removal_reasons(component_index) = append_reason( ...
                removal_reasons(component_index), 'Brain確率未評価');
        elseif ~brain_ok(component_index)
            removal_reasons(component_index) = append_reason( ...
                removal_reasons(component_index), sprintf( ...
                'Brain=%.3f<%.2f', brain_probability(component_index), ...
                config.brain_keep_threshold));
        end
        if ~isfinite(eog_correlations(component_index))
            removal_reasons(component_index) = append_reason( ...
                removal_reasons(component_index), 'EOG相関未評価');
        elseif ~eog_ok(component_index)
            removal_reasons(component_index) = append_reason( ...
                removal_reasons(component_index), sprintf( ...
                'EOG|r|=%.3f>=%.2f', ...
                abs(eog_correlations(component_index)), ...
                config.eog_correlation_threshold));
        end
    end
    if all(remove_mask)
        error('run_pae_specificity_analysis:NoComponentKept', ...
            ['Brain確率50%以上かつ|IC-EOG相関|<0.40を満たすICが' ...
             'ありません。ICA QCを確認してください。']);
    end
end

function output = append_reason(input, addition)
    if strlength(input) == 0
        output = string(addition);
    else
        output = input + "; " + string(addition);
    end
end

function save_joint_ica_qc(merged, training, report, ...
        qc_folder, config)
    n_components = size(merged.icaweights, 1);
    component = (1:n_components)';
    dominant_class = strings(n_components, 1);
    dominant_probability = nan(n_components, 1);
    for component_index = 1:n_components
        probabilities = report.Classification(component_index, :);
        if all(~isfinite(probabilities))
            dominant_class(component_index) = "未評価";
        else
            [dominant_probability(component_index), maximum_index] = ...
                max(probabilities, [], 'omitnan');
            dominant_class(component_index) = ...
                report.Classes(maximum_index);
        end
    end

    kept = ismember(component, report.KeptComponents(:));
    review = ismember(component, report.ReviewComponents(:));
    removed = ismember(component, report.RemovedComponents(:));
    report_table = table(component, dominant_class, ...
        dominant_probability, report.EOGCorrelations(:), kept, review, ...
        removed, report.ReviewReason(:), report.RemovalReason(:), ...
        'VariableNames', {'IC', 'DominantClass', ...
        'DominantProbability', 'EOGCorrelation', 'Kept', 'Review', ...
        'Removed', 'ReviewReason', 'RemovalReason'});
    for class_index = 1:numel(report.Classes)
        variable_name = matlab.lang.makeUniqueStrings( ...
            matlab.lang.makeValidName(report.Classes(class_index)), ...
            report_table.Properties.VariableNames);
        report_table.(variable_name) = ...
            report.Classification(:, class_index);
    end
    writetable(report_table, ...
        fullfile(qc_folder, '共通ICA_ICLabel_EOG判定.csv'), ...
        'Encoding', 'UTF-8');

    icaweights = merged.icaweights;
    icasphere = merged.icasphere;
    icawinv = merged.icawinv;
    icachansind = merged.icachansind;
    channel_labels = string({merged.chanlocs.labels});
    analysis_config = config;
    save(fullfile(qc_folder, '共通ICA_分解と除去判定.mat'), ...
        'icaweights', 'icasphere', 'icawinv', 'icachansind', ...
        'channel_labels', 'report', 'analysis_config', '-v7');

    try
        component_activity = double(training.icaweights * ...
            training.icasphere * ...
            training.data(training.icachansind, :));
        [spectrum, spectrum_freqs] = spectopo(component_activity, ...
            0, training.srate, 'plot', 'off');
        figure_height = max(700, 235 * n_components);
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [50 20 1450 figure_height]);
        cleanup = onCleanup(@() close(fig));
        for component_index = 1:n_components
            ax_map = subplot(n_components, 2, ...
                2 * component_index - 1, 'Parent', fig);
            axes(ax_map); %#ok<LAXES>
            topoplot(merged.icawinv(:, component_index), ...
                merged.chanlocs(merged.icachansind), ...
                'electrodes', 'on', 'style', 'both', ...
                'headrad', 0.5);
            decision_mark = " [Brain50採用]";
            if ismember(component_index, report.RemovedComponents)
                decision_mark = " [除去]";
            end
            map_title = title(ax_map, sprintf( ...
                'IC %d: %s %.1f%% / EOG r=%.3f%s', ...
                component_index, dominant_class(component_index), ...
                100 * dominant_probability(component_index), ...
                report.EOGCorrelations(component_index), decision_mark), ...
                'Interpreter', 'none');
            map_title.Color = 'k';
            set(ax_map, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');

            ax_spectrum = subplot(n_components, 2, ...
                2 * component_index, 'Parent', fig);
            plot(ax_spectrum, spectrum_freqs, ...
                spectrum(component_index, :), 'LineWidth', 1.3);
            xlim(ax_spectrum, [1, min(config.ica_lowpass_hz, ...
                training.srate / 2)]);
            xline(ax_spectrum, 50, '--r', '50 Hz');
            grid(ax_spectrum, 'on');
            xlabel(ax_spectrum, '周波数 (Hz)');
            ylabel(ax_spectrum, 'Power (dB)');
            spectrum_title = title(ax_spectrum, ...
                sprintf('IC %d spectrum', ...
                component_index));
            spectrum_title.Color = 'k';
            set(ax_spectrum, 'Color', 'w', ...
                'XColor', 'k', 'YColor', 'k');
        end
        overall_title = sgtitle(fig, sprintf( ...
            ['PAE特異性感度分析: 採用IC = %s / 除去IC = %s' ...
             '（Brain >= %.0f%% かつ |IC-EOG r| < %.2f）'], ...
            mat2str(report.KeptComponents), ...
            mat2str(report.RemovedComponents), ...
            100 * config.brain_keep_threshold, ...
            config.eog_correlation_threshold), ...
            'Interpreter', 'none');
        overall_title.Color = 'k';
        exportgraphics(fig, ...
            fullfile(qc_folder, '共通ICA_ICLabel_EOG判定.png'), ...
            'Resolution', 180);
        clear cleanup;
    catch ME
        warning('run_joint_ica_cwt_band_analysis:QCFigureFailed', ...
            '共通ICA QC画像を保存できませんでした: %s', ME.message);
    end
end

function [manifest, results] = analyze_all_conditions( ...
        datasets, output_root, joint_report, config)
    n_conditions = numel(datasets);
    n_channels = numel(config.analysis_channels);
    results = repmat(struct(), n_conditions, 1);
    manifest = table();

    for condition_index = 1:n_conditions
        dataset = datasets(condition_index);
        EEG = dataset.EEG;
        event_latencies = get_event_latencies(EEG, config.event_type);
        if numel(event_latencies) < 2
            error('run_joint_ica_cwt_band_analysis:TooFewEvents', ...
                '%s のSaccadeイベントが不足しています。', ...
                dataset.SourceSET);
        end
        EEG_epoch = epoch_for_high_frequency(EEG, config);
        eog_metrics = calculate_epoch_eog_metrics( ...
            EEG_epoch, dataset.EOGIndex, config);
        fprintf('  ERSP: %s（低周波候補%d / 高周波epoch%d）\n', ...
            dataset.ConditionLabel, numel(event_latencies), ...
            EEG_epoch.trials);

        results(condition_index).ConditionLabel = ...
            dataset.ConditionLabel;
        results(condition_index).ConditionValue = ...
            dataset.ConditionValue;
        results(condition_index).IsGreen = dataset.IsGreen;
        results(condition_index).SourceSET = dataset.SourceSET;
        results(condition_index).BaseName = dataset.BaseName;
        results(condition_index).EOGMetrics = eog_metrics;
        channel_template = struct('Name', "", 'Mean', [], ...
            'SEM', [], 'CI95Lower', [], 'CI95Upper', [], ...
            'N', [], 'Trials', [], 'InducedMean', [], ...
            'InducedSEM', [], 'InducedN', [], 'InducedTrials', [], ...
            'Evoked', [], 'ITC', []);
        results(condition_index).Channels = repmat(channel_template, ...
            n_channels, 1);

        for channel_position = 1:n_channels
            channel_name = config.analysis_channels(channel_position);
            target_position = find(strcmpi(config.target_channels, ...
                channel_name), 1);
            channel_index = dataset.TargetIndices(target_position);
            fprintf('    CWT/帯域: %s ...\n', channel_name);
            [low_trials, low_freqs, low_valid_count] = ...
                compute_continuous_trial_ersp(EEG, channel_index, ...
                event_latencies, config);
            [high_trials, high_induced_trials, high_evoked, ...
                high_itc, high_freqs] = ...
                compute_epoched_trial_decomposition(EEG_epoch, ...
                    channel_index, config);

            [low_mean, low_sem, low_n] = ...
                mean_sem_over_trials(low_trials);
            [high_mean, high_sem, high_n] = ...
                mean_sem_over_trials(high_trials);
            ersp_mean = [low_mean; high_mean];
            ersp_sem = [low_sem; high_sem];
            ersp_n = [low_n; high_n];
            freqs = [low_freqs(:); high_freqs(:)];
            times = config.common_times_ms(:)';

            band_stats = calculate_trial_band_statistics( ...
                low_trials, low_freqs, high_trials, high_freqs, ...
                times, config);
            induced_band_stats = ...
                calculate_high_frequency_band_statistics( ...
                    high_induced_trials, high_freqs, times, config);
            evoked_band = calculate_high_frequency_band_curves( ...
                high_evoked, high_freqs, config);
            itc_band = calculate_high_frequency_band_curves( ...
                high_itc, high_freqs, config);

            channel_folder = fullfile(output_root, channel_name);
            ensure_folder(channel_folder);
            safe_base = safe_filename(dataset.BaseName);
            file_stem = safe_base + "_" + channel_name;
            cwt_mat = fullfile(channel_folder, ...
                file_stem + "_共通ICA_CWT.mat");
            cwt_png = fullfile(channel_folder, ...
                file_stem + "_共通ICA_CWT.png");
            band_mat = fullfile(channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.mat");
            band_png = fullfile(channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.png");
            band_csv = fullfile(channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.csv");
            decomposition_png = fullfile(channel_folder, ...
                file_stem + "_Brain50_高周波成分分解.png");

            source_set = dataset.SourceSET;
            condition_label = dataset.ConditionLabel;
            condition_value_hz = dataset.ConditionValue;
            is_green_noled = dataset.IsGreen;
            kept_components = joint_report.KeptComponents;
            review_components = joint_report.ReviewComponents;
            removed_components = joint_report.RemovedComponents;
            analysis_config = config;
            low_trial_count = size(low_trials, 3);
            high_trial_count = size(high_trials, 3);
            save(cwt_mat, 'source_set', 'condition_label', ...
                'condition_value_hz', 'is_green_noled', ...
                'kept_components', 'review_components', ...
                'removed_components', ...
                'analysis_config', ...
                'ersp_mean', 'ersp_sem', 'ersp_n', 'freqs', ...
                'times', 'low_freqs', 'high_freqs', ...
                'low_trial_count', 'high_trial_count', '-v7');

            band_mean = band_stats.Mean;
            band_sem = band_stats.SEM;
            band_ci95_lower = band_stats.CI95Lower;
            band_ci95_upper = band_stats.CI95Upper;
            band_n = band_stats.N;
            band_trials = band_stats.Trials;
            band_summary = band_stats.Summary;
            induced_band_mean = induced_band_stats.Mean;
            induced_band_sem = induced_band_stats.SEM;
            induced_band_n = induced_band_stats.N;
            induced_band_trials = induced_band_stats.Trials;
            evoked_band_ersp = evoked_band;
            band_itc = itc_band;
            save(band_mat, 'source_set', 'condition_label', ...
                'condition_value_hz', 'is_green_noled', ...
                'kept_components', 'review_components', ...
                'removed_components', ...
                'analysis_config', 'times', ...
                'band_mean', 'band_sem', 'band_ci95_lower', ...
                'band_ci95_upper', 'band_n', 'band_trials', ...
                'band_summary', 'induced_band_mean', ...
                'induced_band_sem', 'induced_band_n', ...
                'induced_band_trials', 'evoked_band_ersp', ...
                'band_itc', 'eog_metrics', '-v7');
            band_table = build_band_timecourse_table( ...
                band_stats, times, dataset, channel_name, config);
            writetable(band_table, band_csv, 'Encoding', 'UTF-8');

            create_cwt_mean_figure(times, freqs, ersp_mean, ...
                dataset.ConditionLabel, channel_name, cwt_png, config);
            create_single_condition_band_figure(times, band_stats, ...
                dataset.ConditionLabel, channel_name, band_png, config);
            create_decomposition_figure(times, band_stats.Mean, ...
                induced_band_stats.Mean, evoked_band, itc_band, ...
                dataset.ConditionLabel, channel_name, ...
                decomposition_png, config);

            channel_result = struct();
            channel_result.Name = channel_name;
            channel_result.Mean = band_stats.Mean;
            channel_result.SEM = band_stats.SEM;
            channel_result.CI95Lower = band_stats.CI95Lower;
            channel_result.CI95Upper = band_stats.CI95Upper;
            channel_result.N = band_stats.N;
            channel_result.Trials = band_stats.Trials;
            channel_result.InducedMean = induced_band_stats.Mean;
            channel_result.InducedSEM = induced_band_stats.SEM;
            channel_result.InducedN = induced_band_stats.N;
            channel_result.InducedTrials = induced_band_stats.Trials;
            channel_result.Evoked = evoked_band;
            channel_result.ITC = itc_band;
            results(condition_index).Channels(channel_position) = ...
                channel_result;

            row = table(dataset.SourceSET, dataset.ConditionLabel, ...
                dataset.ConditionValue, dataset.IsGreen, channel_name, ...
                dataset.RawSaccadeEvents, dataset.DuplicatesRemoved, ...
                low_valid_count, EEG_epoch.trials, ...
                string(mat2str(joint_report.KeptComponents)), ...
                string(mat2str(joint_report.ReviewComponents)), ...
                string(mat2str(joint_report.RemovedComponents)), ...
                string(cwt_mat), string(cwt_png), string(band_mat), ...
                string(band_png), string(band_csv), ...
                string(decomposition_png), ...
                'VariableNames', {'SourceSET', 'Condition', ...
                'ConditionValueHz', 'IsGreenNoLED', 'Channel', ...
                'RawSaccadeEvents', 'DuplicateEventsRemoved', ...
                'LowFrequencyValidTrials', 'HighFrequencyTrials', ...
                'CommonKeptICs', 'CommonReviewICs', ...
                'CommonRemovedICs', 'CWTMAT', 'CWTPNG', ...
                'BandMAT', 'BandPNG', 'BandCSV', ...
                'HighFrequencyDecompositionPNG'});
            manifest = append_table_row(manifest, row);

            clear low_trials high_trials high_induced_trials ...
                band_trials induced_band_trials
        end
    end
end

function latencies = get_event_latencies(EEG, event_type)
    types = event_types_as_strings(EEG.event);
    event_indices = find(strcmpi(types, event_type));
    latencies = double([EEG.event(event_indices).latency]);
    latencies = latencies(isfinite(latencies));
end

function EEG_epoch = epoch_for_high_frequency(EEG, config)
    EEG_epoch = pop_epoch(EEG, {char(config.event_type)}, ...
        config.high_epoch_window_s, 'epochinfo', 'yes');
    if EEG_epoch.trials < 2
        error('run_joint_ica_cwt_band_analysis:TooFewHighEpochs', ...
            '端部を除外すると高周波ERSP用epochが2試行未満です。');
    end
    EEG_epoch = eeg_checkset(EEG_epoch, 'eventconsistency');
end

function [trial_ersp, freqs, valid_event_count] = ...
        compute_continuous_trial_ersp(EEG, channel_index, ...
        event_latencies, config)
    signal = double(reshape(EEG.data(channel_index, :, :), [], 1));
    requested_freqs = config.low_frequencies_hz;
    longest_wavelet_samples = ceil( ...
        config.wavelet_cycles * EEG.srate / min(requested_freqs));
    requested_time_points = floor( ...
        (EEG.xmax - EEG.xmin) * 1000 / config.output_step_ms) + 1;
    maximum_time_points = EEG.pnts - longest_wavelet_samples;
    n_timesout = min(requested_time_points, maximum_time_points);
    if n_timesout < 10
        error('run_joint_ica_cwt_band_analysis:RecordingTooShort', ...
            '0.5 Hz・3サイクルCWTに対して連続記録が短すぎます。');
    end

    [tf, freqs, continuous_times_ms] = timefreq( ...
        signal, EEG.srate, 'cycles', config.wavelet_cycles, ...
        'freqs', requested_freqs, 'ntimesout', n_timesout, ...
        'tlimits', [EEG.xmin EEG.xmax] * 1000, ...
        'padratio', config.padratio, 'verbose', 'off');
    tf = squeeze(tf);
    if size(tf, 1) ~= numel(freqs)
        tf = tf.';
    end
    log_power = 10 * log10(max(abs(tf).^2, realmin('double')));

    display_times = config.common_times_ms;
    baseline_times = config.baseline_window_ms(1):...
        config.output_step_ms:config.baseline_window_ms(2);
    n_events = numel(event_latencies);
    trial_ersp = nan(numel(freqs), numel(display_times), n_events);
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
        display_log_power = interp1(continuous_times_ms(:), ...
            log_power.', display_query(:), 'linear').';
        baseline_log_power = interp1(continuous_times_ms(:), ...
            log_power.', baseline_query(:), 'linear').';
        baseline_count = sum(isfinite(baseline_log_power), 2);
        baseline_reference = sum(baseline_log_power, 2, ...
            'omitnan') ./ max(baseline_count, 1);
        trial_ersp(:, :, event_index) = ...
            display_log_power - baseline_reference;
        valid(event_index) = true;
    end

    trial_ersp = trial_ersp(:, :, valid);
    valid_event_count = size(trial_ersp, 3);
    if valid_event_count < 2
        error('run_joint_ica_cwt_band_analysis:TooFewLowEvents', ...
            '低周波CWTの端部条件を満たすSaccadeが2試行未満です。');
    end
    freqs = double(freqs(:));
end

function metrics = calculate_epoch_eog_metrics( ...
        EEG_epoch, eog_index, config)
    eog = double(reshape(EEG_epoch.data(eog_index, :, :), ...
        EEG_epoch.pnts, EEG_epoch.trials));
    times_ms = double(EEG_epoch.times(:));
    metric_mask = times_ms >= config.eog_metric_window_ms(1) & ...
        times_ms <= config.eog_metric_window_ms(2);
    if sum(metric_mask) < 2
        error('run_pae_specificity_analysis:MissingEOGMetricWindow', ...
            'EOG試行指標の時間窓がepoch内にありません。');
    end
    windowed = eog(metric_mask, :);
    peak_to_peak = max(windowed, [], 1, 'omitnan') - ...
        min(windowed, [], 1, 'omitnan');
    rms_value = sqrt(mean(windowed .^ 2, 1, 'omitnan'));
    derivative = diff(windowed, 1, 1) * EEG_epoch.srate;
    derivative_peak = max(abs(derivative), [], 1, 'omitnan');
    metrics = table((1:EEG_epoch.trials)', peak_to_peak(:), ...
        rms_value(:), derivative_peak(:), ...
        'VariableNames', {'Trial', 'EOGPeakToPeak', 'EOGRMS', ...
        'EOGDerivativePeak'});
end

function [trial_ersp_common, induced_ersp_common, ...
        evoked_ersp_common, itc_common, freqs] = ...
        compute_epoched_trial_decomposition( ...
        EEG_epoch, channel_index, config)
    signal = double(reshape(EEG_epoch.data(channel_index, :, :), ...
        EEG_epoch.pnts, EEG_epoch.trials));
    [tf, freqs, times_ms] = timefreq(signal, EEG_epoch.srate, ...
        'cycles', config.wavelet_cycles, ...
        'freqs', config.high_frequencies_hz, ...
        'ntimesout', config.high_timesout, ...
        'tlimits', [EEG_epoch.xmin EEG_epoch.xmax] * 1000, ...
        'padratio', config.padratio, 'verbose', 'off');
    tf = normalize_tf_dimensions(tf, numel(freqs), ...
        EEG_epoch.trials);
    trial_ersp = baseline_normalize_log_power(tf, times_ms, config);

    phase_unit = tf ./ max(abs(tf), realmin('double'));
    itc = abs(mean(phase_unit, 3, 'omitnan'));

    evoked_signal = mean(signal, 2, 'omitnan');
    induced_signal = signal - evoked_signal;
    [induced_tf, induced_freqs, induced_times_ms] = timefreq( ...
        induced_signal, EEG_epoch.srate, ...
        'cycles', config.wavelet_cycles, ...
        'freqs', config.high_frequencies_hz, ...
        'ntimesout', config.high_timesout, ...
        'tlimits', [EEG_epoch.xmin EEG_epoch.xmax] * 1000, ...
        'padratio', config.padratio, 'verbose', 'off');
    induced_tf = normalize_tf_dimensions(induced_tf, ...
        numel(induced_freqs), EEG_epoch.trials);
    induced_trial_ersp = baseline_normalize_log_power( ...
        induced_tf, induced_times_ms, config);

    [evoked_tf, evoked_freqs, evoked_times_ms] = timefreq( ...
        evoked_signal, EEG_epoch.srate, ...
        'cycles', config.wavelet_cycles, ...
        'freqs', config.high_frequencies_hz, ...
        'ntimesout', config.high_timesout, ...
        'tlimits', [EEG_epoch.xmin EEG_epoch.xmax] * 1000, ...
        'padratio', config.padratio, 'verbose', 'off');
    evoked_tf = normalize_tf_dimensions(evoked_tf, ...
        numel(evoked_freqs), 1);
    evoked_trial_ersp = baseline_normalize_log_power( ...
        evoked_tf, evoked_times_ms, config);

    if numel(induced_freqs) ~= numel(freqs) || ...
            numel(evoked_freqs) ~= numel(freqs) || ...
            any(abs(double(induced_freqs(:)) - double(freqs(:))) > 1e-9) || ...
            any(abs(double(evoked_freqs(:)) - double(freqs(:))) > 1e-9)
        error('run_pae_specificity_analysis:FrequencyMismatch', ...
            'Total/Induced/Evokedの周波数グリッドが一致しません。');
    end

    trial_ersp_common = interpolate_tf_trials(trial_ersp, ...
        times_ms, config.common_times_ms);
    induced_ersp_common = interpolate_tf_trials( ...
        induced_trial_ersp, induced_times_ms, ...
        config.common_times_ms);
    evoked_interpolated = interpolate_tf_trials( ...
        evoked_trial_ersp, evoked_times_ms, ...
        config.common_times_ms);
    evoked_ersp_common = evoked_interpolated(:, :, 1);
    itc_interpolated = interpolate_tf_trials(reshape(itc, ...
        size(itc, 1), size(itc, 2), 1), times_ms, ...
        config.common_times_ms);
    itc_common = itc_interpolated(:, :, 1);
    freqs = double(freqs(:));
end

function trial_ersp = baseline_normalize_log_power( ...
        tf, times_ms, config)
    log_power = 10 * log10(max(abs(tf).^2, realmin('double')));
    baseline_mask = times_ms >= config.baseline_window_ms(1) & ...
        times_ms <= config.baseline_window_ms(2);
    if sum(baseline_mask) < 2
        error('run_joint_ica_cwt_band_analysis:MissingBaseline', ...
            '高周波ERSPのベースライン点が不足しています。');
    end
    baseline_reference = mean(log_power(:, baseline_mask, :), ...
        2, 'omitnan');
    trial_ersp = log_power - baseline_reference;
end

function trial_common = interpolate_tf_trials( ...
        trial_values, times_ms, common_times_ms)
    n_frequencies = size(trial_values, 1);
    n_trials = size(trial_values, 3);
    trial_common = nan(n_frequencies, numel(common_times_ms), n_trials);
    for trial_index = 1:n_trials
        interpolated = interp1(double(times_ms(:)), ...
            double(trial_values(:, :, trial_index).'), ...
            common_times_ms(:), 'linear', NaN);
        trial_common(:, :, trial_index) = interpolated.';
    end
end

function tf = normalize_tf_dimensions(tf, n_frequencies, n_trials)
    tf_size = size(tf);
    if n_trials == 1 && ismatrix(tf)
        if size(tf, 1) ~= n_frequencies
            tf = tf.';
        end
        tf = reshape(tf, n_frequencies, size(tf, 2), 1);
        return;
    end
    if size(tf, 1) == n_frequencies && size(tf, 3) == n_trials
        return;
    end
    if size(tf, 2) == n_frequencies && size(tf, 3) == n_trials
        tf = permute(tf, [2 1 3]);
        return;
    end
    error('run_joint_ica_cwt_band_analysis:UnexpectedTFSize', ...
        'timefreq出力サイズが想定外です: %s', mat2str(tf_size));
end

function [mean_values, sem_values, n_values] = ...
        mean_sem_over_trials(trial_values)
    mean_values = mean(trial_values, 3, 'omitnan');
    n_values = sum(isfinite(trial_values), 3);
    standard_deviation = std(trial_values, 0, 3, 'omitnan');
    sem_values = standard_deviation ./ sqrt(max(n_values, 1));
    sem_values(n_values < 2) = NaN;
end

function stats = calculate_trial_band_statistics( ...
        low_trials, low_freqs, high_trials, high_freqs, times, config)
    n_bands = numel(config.band_names);
    n_times = numel(times);
    stats.Mean = nan(n_bands, n_times);
    stats.SEM = nan(n_bands, n_times);
    stats.CI95Lower = nan(n_bands, n_times);
    stats.CI95Upper = nan(n_bands, n_times);
    stats.N = zeros(n_bands, n_times);
    stats.Trials = cell(n_bands, 1);
    coverage = strings(n_bands, 1);
    frequency_bin_count = zeros(n_bands, 1);
    excluded_line_bin_count = zeros(n_bands, 1);
    actual_minimum_hz = nan(n_bands, 1);
    actual_maximum_hz = nan(n_bands, 1);

    for band_index = 1:n_bands
        low_hz = config.band_ranges(band_index, 1);
        high_hz = config.band_ranges(band_index, 2);
        if high_hz <= max(low_freqs)
            source_freqs = low_freqs;
            source_trials = low_trials;
        else
            source_freqs = high_freqs;
            source_trials = high_trials;
        end
        band_mask = source_freqs >= low_hz & source_freqs <= high_hz;
        excluded_line_bin_count(band_index) = 0;
        frequency_bin_count(band_index) = sum(band_mask);
        if ~any(band_mask)
            coverage(band_index) = "なし";
            stats.Trials{band_index} = nan(n_times, 0);
            continue;
        end

        actual_minimum_hz(band_index) = min(source_freqs(band_mask));
        actual_maximum_hz(band_index) = max(source_freqs(band_mask));
        if min(source_freqs) <= low_hz && ...
                max(source_freqs) >= high_hz
            coverage(band_index) = "全域";
        else
            coverage(band_index) = "一部";
        end
        values = mean(source_trials(band_mask, :, :), 1, 'omitnan');
        values = reshape(values, n_times, size(source_trials, 3));
        stats.Trials{band_index} = values;
        stats.Mean(band_index, :) = mean(values, 2, 'omitnan').';
        n_valid = sum(isfinite(values), 2);
        stats.N(band_index, :) = n_valid.';
        standard_deviation = std(values, 0, 2, 'omitnan');
        sem = standard_deviation ./ sqrt(max(n_valid, 1));
        sem(n_valid < 2) = NaN;
        stats.SEM(band_index, :) = sem.';
        stats.CI95Lower(band_index, :) = ...
            stats.Mean(band_index, :) - ...
            config.confidence_multiplier * stats.SEM(band_index, :);
        stats.CI95Upper(band_index, :) = ...
            stats.Mean(band_index, :) + ...
            config.confidence_multiplier * stats.SEM(band_index, :);
    end

    stats.Summary = table(config.band_names, config.band_labels, ...
        config.band_ranges(:, 1), config.band_ranges(:, 2), ...
        coverage, config.band_interpretation, frequency_bin_count, ...
        excluded_line_bin_count, actual_minimum_hz, ...
        actual_maximum_hz, ...
        'VariableNames', {'Band', 'BandLabelJapanese', ...
        'DefinedMinimumHz', 'DefinedMaximumHz', 'Coverage', ...
        'Interpretation', 'FrequencyBinCount', ...
        'ExcludedLineNoiseBinCount', 'ActualMinimumHz', ...
        'ActualMaximumHz'});
end

function stats = calculate_high_frequency_band_statistics( ...
        high_trials, high_freqs, times, config)
    n_bands = numel(config.band_names);
    n_times = numel(times);
    stats.Mean = nan(n_bands, n_times);
    stats.SEM = nan(n_bands, n_times);
    stats.CI95Lower = nan(n_bands, n_times);
    stats.CI95Upper = nan(n_bands, n_times);
    stats.N = zeros(n_bands, n_times);
    stats.Trials = cell(n_bands, 1);
    for band_index = 1:n_bands
        low_hz = config.band_ranges(band_index, 1);
        high_hz = config.band_ranges(band_index, 2);
        band_mask = high_freqs >= low_hz & high_freqs <= high_hz;
        if ~any(band_mask)
            stats.Trials{band_index} = nan(n_times, 0);
            continue;
        end
        values = mean(high_trials(band_mask, :, :), 1, 'omitnan');
        values = reshape(values, n_times, size(high_trials, 3));
        stats.Trials{band_index} = values;
        stats.Mean(band_index, :) = mean(values, 2, 'omitnan').';
        n_valid = sum(isfinite(values), 2);
        stats.N(band_index, :) = n_valid.';
        standard_deviation = std(values, 0, 2, 'omitnan');
        sem = standard_deviation ./ sqrt(max(n_valid, 1));
        sem(n_valid < 2) = NaN;
        stats.SEM(band_index, :) = sem.';
        stats.CI95Lower(band_index, :) = ...
            stats.Mean(band_index, :) - ...
            config.confidence_multiplier * stats.SEM(band_index, :);
        stats.CI95Upper(band_index, :) = ...
            stats.Mean(band_index, :) + ...
            config.confidence_multiplier * stats.SEM(band_index, :);
    end
end

function band_curves = calculate_high_frequency_band_curves( ...
        frequency_time_values, high_freqs, config)
    band_curves = nan(numel(config.band_names), ...
        size(frequency_time_values, 2));
    for band_index = 1:numel(config.band_names)
        band_mask = high_freqs >= config.band_ranges(band_index, 1) & ...
            high_freqs <= config.band_ranges(band_index, 2);
        if any(band_mask)
            band_curves(band_index, :) = mean( ...
                frequency_time_values(band_mask, :), 1, 'omitnan');
        end
    end
end

function output = append_table_row(output, row)
    if isempty(output)
        output = row;
    else
        output = [output; row];
    end
end

function safe_name = safe_filename(input_name)
    safe_name = regexprep(string(input_name), '[<>:"/\\|?*]', '_');
    safe_name = strip(safe_name);
    if strlength(safe_name) == 0
        safe_name = "unnamed";
    end
end

function output_table = build_band_timecourse_table( ...
        stats, times, dataset, channel_name, config)
    n_bands = numel(config.band_names);
    n_times = numel(times);
    source_set = repmat(dataset.SourceSET, n_bands * n_times, 1);
    condition = repmat(dataset.ConditionLabel, ...
        n_bands * n_times, 1);
    channel = repmat(string(channel_name), n_bands * n_times, 1);
    band = repelem(config.band_names, n_times);
    band_label = repelem(config.band_labels, n_times);
    interpretation = repelem(config.band_interpretation, n_times);
    time_ms = repmat(times(:), n_bands, 1);
    mean_ersp_db = reshape(stats.Mean.', [], 1);
    sem_db = reshape(stats.SEM.', [], 1);
    ci95_lower_db = reshape(stats.CI95Lower.', [], 1);
    ci95_upper_db = reshape(stats.CI95Upper.', [], 1);
    trial_count = reshape(stats.N.', [], 1);
    output_table = table(source_set, condition, channel, band, ...
        band_label, interpretation, time_ms, mean_ersp_db, sem_db, ...
        ci95_lower_db, ci95_upper_db, trial_count, ...
        'VariableNames', {'SourceSET', 'Condition', 'Channel', ...
        'Band', 'BandLabelJapanese', 'Interpretation', 'TimeMs', ...
        'MeanERSPdB', 'TrialSEMDb', 'CI95LowerDb', ...
        'CI95UpperDb', 'TrialCount'});
end

function create_cwt_mean_figure(times, freqs, ersp_mean, ...
        condition_label, channel_name, output_png, config)
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [80 80 1250 720]);
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    surface(ax, times, freqs, zeros(size(ersp_mean)), ersp_mean, ...
        'EdgeColor', 'none', 'FaceColor', 'flat');
    view(ax, 2);
    axis(ax, 'tight');
    xlim(ax, config.display_window_ms);
    ylim(ax, [min(freqs), max(freqs)]);
    colormap(ax, turbo);
    colorbar_handle = colorbar(ax);
    colorbar_handle.Label.String = '試行平均ERSP (dB)';
    colorbar_handle.Color = 'k';
    finite_values = ersp_mean(isfinite(ersp_mean));
    if ~isempty(finite_values)
        symmetric_limit = finite_percentile(abs(finite_values), 98);
        if isfinite(symmetric_limit) && symmetric_limit > 0
            clim(ax, [-symmetric_limit symmetric_limit]);
        end
    end
    hold(ax, 'on');
    xline(ax, 0, '--m', 'LineWidth', 1.4);
    yline(ax, 7.5, ':k', '低周波/高周波CWT境界', ...
        'LabelHorizontalAlignment', 'left');
    xlabel(ax, '時間 (ms)');
    ylabel(ax, '周波数 (Hz)');
    title_handle = title(ax, sprintf( ...
        '%s / %s: 全条件共通ICA後の試行平均ERSP（0 ms = サッケード開始）', ...
        condition_label, channel_name), 'Interpreter', 'none');
    title_handle.Color = 'k';
    set(ax, 'FontSize', 11, 'Box', 'on', ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function value = finite_percentile(values, percentile)
    values = sort(double(values(isfinite(values))));
    if isempty(values)
        value = NaN;
        return;
    end
    position = 1 + (numel(values) - 1) * percentile / 100;
    lower_index = floor(position);
    upper_index = ceil(position);
    if lower_index == upper_index
        value = values(lower_index);
    else
        fraction = position - lower_index;
        value = values(lower_index) * (1 - fraction) + ...
            values(upper_index) * fraction;
    end
end

function create_single_condition_band_figure(times, stats, ...
        condition_label, channel_name, output_png, config)
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [60 30 1400 1220]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 4, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    line_color = [0.05 0.35 0.72];

    for band_index = 1:numel(config.band_names)
        ax = nexttile(layout);
        mean_values = stats.Mean(band_index, :);
        plot_mean_line(ax, times, mean_values, line_color, "");
        hold(ax, 'on');
        xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1);
        yline(ax, 0, ':k', 'LineWidth', 0.8);
        xlim(ax, config.display_window_ms);
        ylim(ax, curve_limits(mean_values, zeros(size(mean_values))));
        grid(ax, 'on');
        xlabel(ax, '時間 (ms)');
        ylabel(ax, '平均ERSP (dB)');
        title_handle = title(ax, sprintf('%s (%.1f–%.1f Hz) [%s]', ...
            config.band_labels(band_index), ...
            config.band_ranges(band_index, 1), ...
            config.band_ranges(band_index, 2), ...
            config.band_interpretation(band_index)), ...
            'Interpreter', 'none');
        title_handle.Color = 'k';
        set(ax, 'FontSize', 10, 'Box', 'on', ...
            'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    end
    layout_title = title(layout, sprintf( ...
        ['%s / %s: 帯域平均ERSP（線 = 試行平均、' ...
        '0 ms = サッケード開始）'], ...
        condition_label, channel_name), 'Interpreter', 'none', ...
        'FontWeight', 'bold');
    layout_title.Color = 'k';
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function create_decomposition_figure(times, total_mean, ...
        induced_mean, evoked_ersp, itc, condition_label, ...
        channel_name, output_png, config)
    band_indices = config.specificity_band_indices;
    n_rows = numel(band_indices);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [40 20 1550 270 * n_rows]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, n_rows, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    total_color = [0.05 0.35 0.72];
    induced_color = [0.85 0.25 0.12];
    evoked_color = [0.20 0.60 0.25];

    for row_index = 1:n_rows
        band_index = band_indices(row_index);
        ax_power = nexttile(layout);
        plot_mean_line(ax_power, times, total_mean(band_index, :), ...
            total_color, "Total");
        plot_mean_line(ax_power, times, induced_mean(band_index, :), ...
            induced_color, "Induced");
        plot_mean_line(ax_power, times, evoked_ersp(band_index, :), ...
            evoked_color, "Evoked");
        xline(ax_power, 0, '--k', 'HandleVisibility', 'off');
        yline(ax_power, 0, ':k', 'HandleVisibility', 'off');
        xlim(ax_power, config.display_window_ms);
        grid(ax_power, 'on');
        xlabel(ax_power, '時間 (ms)');
        ylabel(ax_power, 'ERSP (dB)');
        title(ax_power, sprintf('%s (%.1f–%.1f Hz): power', ...
            config.band_labels(band_index), ...
            config.band_ranges(band_index, 1), ...
            config.band_ranges(band_index, 2)), ...
            'Interpreter', 'none', 'Color', 'k');
        legend(ax_power, 'Location', 'best', 'Interpreter', 'none');
        set(ax_power, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
            'Box', 'on');

        ax_itc = nexttile(layout);
        plot(ax_itc, times, itc(band_index, :), ...
            'LineWidth', 1.7, 'Color', [0.55 0.12 0.65]);
        xline(ax_itc, 0, '--k');
        xlim(ax_itc, config.display_window_ms);
        ylim(ax_itc, [0 1]);
        grid(ax_itc, 'on');
        xlabel(ax_itc, '時間 (ms)');
        ylabel(ax_itc, 'ITC');
        title(ax_itc, sprintf('%s: 位相同期', ...
            config.band_labels(band_index)), ...
            'Interpreter', 'none', 'Color', 'k');
        set(ax_itc, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
            'Box', 'on');
    end
    title(layout, sprintf( ...
        ['%s / %s: Brain50共通ICA後の高周波分解 ' ...
         '（0 ms = サッケード開始）'], ...
        condition_label, channel_name), ...
        'Interpreter', 'none', 'FontWeight', 'bold', 'Color', 'k');
    exportgraphics(fig, output_png, 'Resolution', 180);
    clear cleanup;
end

function plot_handle = plot_mean_line(ax, times, mean_values, ...
        line_color, display_name)
    times = double(times(:)');
    mean_values = double(mean_values(:)');
    hold(ax, 'on');
    if strlength(display_name) == 0
        plot_handle = plot(ax, times, mean_values, ...
            'LineWidth', 1.7, 'Color', line_color, ...
            'HandleVisibility', 'off');
    else
        plot_handle = plot(ax, times, mean_values, ...
            'LineWidth', 1.7, 'Color', line_color, ...
            'DisplayName', display_name);
    end
end

function limits = curve_limits(mean_values, error_values)
    values = [mean_values(:) - error_values(:); ...
        mean_values(:) + error_values(:); 0];
    values = values(isfinite(values));
    if isempty(values)
        limits = [-1 1];
        return;
    end
    limits = [min(values), max(values)];
    if diff(limits) < 1e-9
        limits = limits + [-1 1];
    else
        margin = max(0.1, 0.08 * diff(limits));
        limits = limits + [-margin margin];
    end
end

function create_all_channel_overlays(results, output_root, config)
    if isempty(results)
        return;
    end
    n_conditions = numel(results);
    n_channels = numel(config.analysis_channels);
    condition_labels = strings(n_conditions, 1);
    condition_values_hz = nan(n_conditions, 1);
    is_green_noled = false(n_conditions, 1);
    source_sets = strings(n_conditions, 1);
    for condition_index = 1:n_conditions
        condition_labels(condition_index) = ...
            results(condition_index).ConditionLabel;
        condition_values_hz(condition_index) = ...
            results(condition_index).ConditionValue;
        is_green_noled(condition_index) = ...
            results(condition_index).IsGreen;
        source_sets(condition_index) = results(condition_index).SourceSET;
    end

    for channel_position = 1:n_channels
        channel_name = config.analysis_channels(channel_position);
        channel_folder = fullfile(output_root, channel_name);
        ensure_folder(channel_folder);
        [raw_mean, raw_sem, raw_n] = collect_channel_arrays( ...
            results, channel_position, config);
        raw_png = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.png');
        raw_mat = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.mat');
        raw_csv = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.csv');
        create_overlay_figure(config.common_times_ms, raw_mean, ...
            condition_labels, channel_name, raw_png, ...
            config, false);
        raw_table = build_overlay_table(config.common_times_ms, ...
            raw_mean, raw_sem, raw_n, condition_labels, ...
            source_sets, channel_name, config);
        writetable(raw_table, raw_csv, 'Encoding', 'UTF-8');
        times = config.common_times_ms;
        analysis_config = config;
        save(raw_mat, 'times', 'raw_mean', 'raw_sem', 'raw_n', ...
            'condition_labels', 'condition_values_hz', ...
            'is_green_noled', 'source_sets', 'analysis_config', '-v7');

        green_indices = find(is_green_noled);
        if isempty(green_indices)
            fprintf('    %s: 緑なしがないため差分重ね描きを省略。\n', ...
                channel_name);
            continue;
        elseif numel(green_indices) > 1
            warning('run_joint_ica_cwt_band_analysis:MultipleGreen', ...
                '%sには緑なしが複数あるため差分を一意に作れません。', ...
                output_root);
            continue;
        end

        green_index = green_indices(1);
        condition_indices = find(~is_green_noled);
        n_differences = numel(condition_indices);
        difference_mean = nan(n_differences, ...
            numel(config.band_names), numel(config.common_times_ms));
        difference_sem = difference_mean;
        condition_n = difference_mean;
        green_n = difference_mean;
        difference_labels = strings(n_differences, 1);
        difference_source_sets = strings(n_differences, 1);
        for difference_index = 1:n_differences
            condition_index = condition_indices(difference_index);
            difference_mean(difference_index, :, :) = ...
                raw_mean(condition_index, :, :) - ...
                raw_mean(green_index, :, :);
            difference_sem(difference_index, :, :) = sqrt( ...
                raw_sem(condition_index, :, :) .^ 2 + ...
                raw_sem(green_index, :, :) .^ 2);
            condition_n(difference_index, :, :) = ...
                raw_n(condition_index, :, :);
            green_n(difference_index, :, :) = ...
                raw_n(green_index, :, :);
            difference_labels(difference_index) = ...
                condition_labels(condition_index) + " - 緑なし";
            difference_source_sets(difference_index) = ...
                source_sets(condition_index);
        end
        difference_png = fullfile(channel_folder, ...
            '条件-緑なし_帯域ERSP差分重ね描き_SEM.png');
        difference_mat = fullfile(channel_folder, ...
            '条件-緑なし_帯域ERSP差分重ね描き_SEM.mat');
        difference_csv = fullfile(channel_folder, ...
            '条件-緑なし_帯域ERSP差分重ね描き_SEM.csv');
        create_overlay_figure(config.common_times_ms, ...
            difference_mean, difference_labels, ...
            channel_name, difference_png, config, true);
        difference_table = build_difference_overlay_table( ...
            config.common_times_ms, difference_mean, difference_sem, ...
            condition_n, green_n, difference_labels, ...
            difference_source_sets, source_sets(green_index), ...
            channel_name, config);
        writetable(difference_table, difference_csv, ...
            'Encoding', 'UTF-8');
        green_source_set = source_sets(green_index);
        save(difference_mat, 'times', 'difference_mean', ...
            'difference_sem', 'condition_n', 'green_n', ...
            'difference_labels', 'difference_source_sets', ...
            'green_source_set', 'analysis_config', '-v7');
    end
end

function [mean_array, sem_array, n_array] = ...
        collect_channel_arrays(results, channel_position, config)
    n_conditions = numel(results);
    array_size = [n_conditions, numel(config.band_names), ...
        numel(config.common_times_ms)];
    mean_array = nan(array_size);
    sem_array = nan(array_size);
    n_array = zeros(array_size);
    for condition_index = 1:n_conditions
        channel = results(condition_index).Channels(channel_position);
        mean_array(condition_index, :, :) = channel.Mean;
        sem_array(condition_index, :, :) = channel.SEM;
        n_array(condition_index, :, :) = channel.N;
    end
end

function create_overlay_figure(times, mean_array, labels, ...
        channel_name, output_png, config, is_difference)
    n_curves = size(mean_array, 1);
    colors = lines(max(n_curves, 4));
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [40 20 1500 1250]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 4, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    for band_index = 1:numel(config.band_names)
        ax = nexttile(layout);
        all_mean = [];
        for curve_index = 1:n_curves
            mean_values = reshape(mean_array(curve_index, ...
                band_index, :), 1, []);
            plot_mean_line(ax, times, mean_values, ...
                colors(curve_index, :), labels(curve_index));
            all_mean = [all_mean; mean_values]; %#ok<AGROW>
        end
        hold(ax, 'on');
        xline(ax, 0, '--', 'Color', [0.25 0.25 0.25], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
        yline(ax, 0, ':k', 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        xlim(ax, config.display_window_ms);
        ylim(ax, curve_limits(all_mean, zeros(size(all_mean))));
        grid(ax, 'on');
        xlabel(ax, '時間 (ms)');
        if is_difference
            ylabel(ax, 'ERSP差分 (dB)');
        else
            ylabel(ax, '平均ERSP (dB)');
        end
        title_handle = title(ax, sprintf('%s (%.1f–%.1f Hz) [%s]', ...
            config.band_labels(band_index), ...
            config.band_ranges(band_index, 1), ...
            config.band_ranges(band_index, 2), ...
            config.band_interpretation(band_index)), ...
            'Interpreter', 'none');
        title_handle.Color = 'k';
        legend_handle = legend(ax, 'Location', 'best', ...
            'Interpreter', 'none', 'FontSize', 8);
        set(legend_handle, 'Color', 'w', 'TextColor', 'k', ...
            'EdgeColor', [0.2 0.2 0.2]);
        set(ax, 'FontSize', 10, 'Box', 'on', ...
            'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    end
    if is_difference
        quantity = '条件 - 緑なし';
    else
        quantity = '0・80・160・緑なし等の全条件';
    end
    layout_title = title(layout, sprintf( ...
        '%s: %s 帯域ERSP重ね描き（線 = 試行平均）', ...
        channel_name, quantity), 'Interpreter', 'none', ...
        'FontWeight', 'bold');
    layout_title.Color = 'k';
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function output_table = build_overlay_table(times, mean_array, ...
        sem_array, n_array, condition_labels, source_sets, ...
        channel_name, config)
    output_table = table();
    for condition_index = 1:size(mean_array, 1)
        for band_index = 1:numel(config.band_names)
            mean_values = reshape(mean_array(condition_index, ...
                band_index, :), [], 1);
            sem_values = reshape(sem_array(condition_index, ...
                band_index, :), [], 1);
            n_values = reshape(n_array(condition_index, ...
                band_index, :), [], 1);
            n_rows = numel(times);
            row = table(repmat(source_sets(condition_index), n_rows, 1), ...
                repmat(condition_labels(condition_index), n_rows, 1), ...
                repmat(string(channel_name), n_rows, 1), ...
                repmat(config.band_names(band_index), n_rows, 1), ...
                repmat(config.band_labels(band_index), n_rows, 1), ...
                times(:), mean_values, sem_values, ...
                mean_values - config.confidence_multiplier * sem_values, ...
                mean_values + config.confidence_multiplier * sem_values, ...
                n_values, ...
                'VariableNames', {'SourceSET', 'Condition', 'Channel', ...
                'Band', 'BandLabelJapanese', 'TimeMs', 'MeanERSPdB', ...
                'TrialSEMDb', 'CI95LowerDb', 'CI95UpperDb', ...
                'TrialCount'});
            output_table = append_table_row(output_table, row);
        end
    end
end

function output_table = build_difference_overlay_table( ...
        times, difference_mean, difference_sem, condition_n, green_n, ...
        labels, condition_sources, green_source, channel_name, config)
    output_table = table();
    for condition_index = 1:size(difference_mean, 1)
        for band_index = 1:numel(config.band_names)
            mean_values = reshape(difference_mean(condition_index, ...
                band_index, :), [], 1);
            sem_values = reshape(difference_sem(condition_index, ...
                band_index, :), [], 1);
            condition_count = reshape(condition_n(condition_index, ...
                band_index, :), [], 1);
            green_count = reshape(green_n(condition_index, ...
                band_index, :), [], 1);
            n_rows = numel(times);
            row = table( ...
                repmat(condition_sources(condition_index), n_rows, 1), ...
                repmat(green_source, n_rows, 1), ...
                repmat(labels(condition_index), n_rows, 1), ...
                repmat(string(channel_name), n_rows, 1), ...
                repmat(config.band_names(band_index), n_rows, 1), ...
                repmat(config.band_labels(band_index), n_rows, 1), ...
                times(:), mean_values, sem_values, ...
                mean_values - config.confidence_multiplier * sem_values, ...
                mean_values + config.confidence_multiplier * sem_values, ...
                condition_count, green_count, ...
                'VariableNames', {'ConditionSourceSET', ...
                'GreenNoLEDSourceSET', 'DifferenceCondition', ...
                'Channel', 'Band', 'BandLabelJapanese', 'TimeMs', ...
                'MeanDifferenceDb', 'IndependentTrialSEMDb', ...
                'CI95LowerDb', 'CI95UpperDb', ...
                'ConditionTrialCount', 'GreenTrialCount'});
            output_table = append_table_row(output_table, row);
        end
    end
end

function create_pae_specificity_outputs(results, output_root, config)
    if isempty(results)
        return;
    end
    n_channels = numel(config.analysis_channels);
    n_bands = numel(config.band_names);
    n_times = numel(config.common_times_ms);
    measure_names = ["Total", "Induced"];
    measure_fields = ["Trials", "InducedTrials"];
    group_labels = ["0 Hz", "80 Hz", "160 Hz", "緑なし"];
    all_trial_table = table();
    all_contrast_table = table();

    for channel_position = 1:n_channels
        channel_name = config.analysis_channels(channel_position);
        channel_folder = fullfile(output_root, channel_name);
        ensure_folder(channel_folder);
        adjusted_mean = nan(numel(measure_names), n_bands, 4, n_times);
        adjusted_lower = adjusted_mean;
        adjusted_upper = adjusted_mean;
        channel_timecourse_table = table();
        channel_contrast_table = table();

        for measure_index = 1:numel(measure_names)
            for band_position = 1:numel(config.specificity_band_indices)
                band_index = config.specificity_band_indices(band_position);
                [trial_curves, group_code, eog_ptp, eog_rms, ...
                    eog_derivative, source_set, condition_label, ...
                    condition_value, is_green, trial_number] = ...
                    gather_specificity_trials(results, channel_position, ...
                        band_index, measure_fields(measure_index));
                if isempty(trial_curves)
                    continue;
                end

                [means, lowers, uppers] = ...
                    fit_eog_adjusted_timecourses(trial_curves, ...
                        group_code, eog_ptp, eog_derivative);
                adjusted_mean(measure_index, band_index, :, :) = ...
                    reshape(means, 1, 1, 4, n_times);
                adjusted_lower(measure_index, band_index, :, :) = ...
                    reshape(lowers, 1, 1, 4, n_times);
                adjusted_upper(measure_index, band_index, :, :) = ...
                    reshape(uppers, 1, 1, 4, n_times);

                timecourse_rows = build_adjusted_timecourse_table( ...
                    means, lowers, uppers, group_labels, channel_name, ...
                    config.band_names(band_index), ...
                    measure_names(measure_index), ...
                    config.common_times_ms);
                channel_timecourse_table = append_table_row( ...
                    channel_timecourse_table, timecourse_rows);

                window_ms = config.specificity_windows_ms( ...
                    band_position, :);
                time_mask = config.common_times_ms >= window_ms(1) & ...
                    config.common_times_ms <= window_ms(2);
                window_value = mean(trial_curves(time_mask, :), ...
                    1, 'omitnan').';
                contrast_rows = fit_eog_adjusted_window_contrasts( ...
                    window_value, group_code, eog_ptp, ...
                    eog_derivative, channel_name, ...
                    config.band_names(band_index), ...
                    measure_names(measure_index), window_ms, config);
                channel_contrast_table = append_table_row( ...
                    channel_contrast_table, contrast_rows);

                n_trials = numel(window_value);
                trial_rows = table(source_set, condition_label, ...
                    condition_value, is_green, ...
                    repmat(channel_name, n_trials, 1), ...
                    repmat(config.band_names(band_index), n_trials, 1), ...
                    repmat(measure_names(measure_index), n_trials, 1), ...
                    repmat(window_ms(1), n_trials, 1), ...
                    repmat(window_ms(2), n_trials, 1), ...
                    trial_number, window_value, eog_ptp, eog_rms, ...
                    eog_derivative, ...
                    'VariableNames', {'SourceSET', 'Condition', ...
                    'ConditionValueHz', 'IsGreenNoLED', 'Channel', ...
                    'Band', 'Method', 'WindowStartMs', 'WindowEndMs', ...
                    'Trial', 'MeanERSPdB', 'EOGPeakToPeak', ...
                    'EOGRMS', 'EOGDerivativePeak'});
                all_trial_table = append_table_row( ...
                    all_trial_table, trial_rows);
            end
        end

        if ~isempty(channel_timecourse_table)
            writetable(channel_timecourse_table, ...
                fullfile(channel_folder, ...
                    'EOG調整_全条件_時系列.csv'), ...
                'Encoding', 'UTF-8');
        end
        if ~isempty(channel_contrast_table)
            writetable(channel_contrast_table, ...
                fullfile(channel_folder, ...
                    'PAE特異性_EOG調整_事前時間窓.csv'), ...
                'Encoding', 'UTF-8');
            all_contrast_table = append_table_row( ...
                all_contrast_table, channel_contrast_table);
        end
        create_eog_adjusted_overlay_figure( ...
            adjusted_mean, group_labels, channel_name, ...
            fullfile(channel_folder, ...
                'EOG調整_全条件_Alpha以上.png'), config, false);
        create_eog_adjusted_overlay_figure( ...
            adjusted_mean, group_labels, channel_name, ...
            fullfile(channel_folder, ...
                'EOG調整_条件-緑なし_Alpha以上.png'), config, true);
        create_specificity_contrast_figure(channel_contrast_table, ...
            channel_name, fullfile(channel_folder, ...
                'PAE特異性_EOG調整_事前時間窓.png'), config);

        analysis_config = config;
        times = config.common_times_ms;
        save(fullfile(channel_folder, ...
            'PAE特異性_EOG調整結果.mat'), 'adjusted_mean', ...
            'adjusted_lower', 'adjusted_upper', 'group_labels', ...
            'times', 'channel_contrast_table', 'analysis_config', '-v7');
    end

    if ~isempty(all_trial_table)
        writetable(all_trial_table, fullfile(output_root, ...
            'PAE特異性_試行別EOGとERSP.csv'), 'Encoding', 'UTF-8');
    end
    if ~isempty(all_contrast_table)
        writetable(all_contrast_table, fullfile(output_root, ...
            'PAE特異性_EOG調整_全チャンネル.csv'), 'Encoding', 'UTF-8');
    end
end

function [trial_curves, group_code, eog_ptp, eog_rms, ...
        eog_derivative, source_set, condition_label, condition_value, ...
        is_green, trial_number] = gather_specificity_trials( ...
        results, channel_position, band_index, measure_field)
    trial_curves = [];
    group_code = [];
    eog_ptp = [];
    eog_rms = [];
    eog_derivative = [];
    source_set = strings(0, 1);
    condition_label = strings(0, 1);
    condition_value = [];
    is_green = false(0, 1);
    trial_number = [];
    for condition_index = 1:numel(results)
        result = results(condition_index);
        code = specificity_group_code(result.ConditionValue, ...
            result.IsGreen);
        if ~isfinite(code)
            continue;
        end
        curves = result.Channels(channel_position).(measure_field) ...
            {band_index};
        metrics = result.EOGMetrics;
        n_trials = min(size(curves, 2), height(metrics));
        if n_trials < 1
            continue;
        end
        trial_curves = [trial_curves, curves(:, 1:n_trials)]; %#ok<AGROW>
        group_code = [group_code; repmat(code, n_trials, 1)]; %#ok<AGROW>
        eog_ptp = [eog_ptp; metrics.EOGPeakToPeak(1:n_trials)]; %#ok<AGROW>
        eog_rms = [eog_rms; metrics.EOGRMS(1:n_trials)]; %#ok<AGROW>
        eog_derivative = [eog_derivative; ...
            metrics.EOGDerivativePeak(1:n_trials)]; %#ok<AGROW>
        source_set = [source_set; ...
            repmat(result.SourceSET, n_trials, 1)]; %#ok<AGROW>
        condition_label = [condition_label; ...
            repmat(result.ConditionLabel, n_trials, 1)]; %#ok<AGROW>
        condition_value = [condition_value; ...
            repmat(result.ConditionValue, n_trials, 1)]; %#ok<AGROW>
        is_green = [is_green; repmat(result.IsGreen, n_trials, 1)]; %#ok<AGROW>
        trial_number = [trial_number; (1:n_trials)']; %#ok<AGROW>
    end
end

function code = specificity_group_code(condition_value, is_green)
    if is_green
        code = 4;
    elseif isfinite(condition_value) && abs(condition_value) < 1e-9
        code = 1;
    elseif isfinite(condition_value) && abs(condition_value - 80) < 1e-9
        code = 2;
    elseif isfinite(condition_value) && abs(condition_value - 160) < 1e-9
        code = 3;
    else
        code = NaN;
    end
end

function [X, valid, standardized_ptp, standardized_derivative] = ...
        build_eog_design(group_code, eog_ptp, eog_derivative)
    valid = isfinite(group_code) & isfinite(eog_ptp) & ...
        isfinite(eog_derivative);
    standardized_ptp = standardize_predictor(eog_ptp);
    standardized_derivative = standardize_predictor(eog_derivative);
    X = [ones(numel(group_code), 1), group_code == 2, ...
        group_code == 3, group_code == 4, standardized_ptp, ...
        standardized_derivative];
end

function standardized = standardize_predictor(values)
    values = double(values(:));
    mean_value = mean(values, 'omitnan');
    standard_deviation = std(values, 0, 'omitnan');
    if ~isfinite(standard_deviation) || standard_deviation < eps
        standardized = zeros(size(values));
        standardized(~isfinite(values)) = NaN;
    else
        standardized = (values - mean_value) / standard_deviation;
    end
end

function [adjusted_mean, ci_lower, ci_upper] = ...
        fit_eog_adjusted_timecourses(trial_curves, group_code, ...
        eog_ptp, eog_derivative)
    n_times = size(trial_curves, 1);
    adjusted_mean = nan(4, n_times);
    ci_lower = adjusted_mean;
    ci_upper = adjusted_mean;
    [X, base_valid] = build_eog_design( ...
        group_code, eog_ptp, eog_derivative);
    prediction_design = [ ...
        1 0 0 0 0 0; ...
        1 1 0 0 0 0; ...
        1 0 1 0 0 0; ...
        1 0 0 1 0 0];
    for time_index = 1:n_times
        y = double(trial_curves(time_index, :).');
        valid = base_valid & isfinite(y);
        [beta, covariance] = fit_ols(X(valid, :), y(valid));
        if isempty(beta)
            continue;
        end
        predicted = prediction_design * beta;
        prediction_variance = diag( ...
            prediction_design * covariance * prediction_design.');
        prediction_se = sqrt(max(prediction_variance, 0));
        adjusted_mean(:, time_index) = predicted;
        ci_lower(:, time_index) = predicted - 1.96 * prediction_se;
        ci_upper(:, time_index) = predicted + 1.96 * prediction_se;
    end
end

function [beta, covariance, degrees_of_freedom] = fit_ols(X, y)
    beta = [];
    covariance = [];
    degrees_of_freedom = NaN;
    if size(X, 1) <= size(X, 2) || rank(X) < 4
        return;
    end
    beta = pinv(X) * y;
    residual = y - X * beta;
    degrees_of_freedom = size(X, 1) - rank(X);
    if degrees_of_freedom < 1
        beta = [];
        return;
    end
    residual_variance = sum(residual .^ 2) / degrees_of_freedom;
    covariance = residual_variance * pinv(X.' * X);
end

function output_table = fit_eog_adjusted_window_contrasts( ...
        y, group_code, eog_ptp, eog_derivative, channel_name, ...
        band_name, measure_name, window_ms, config)
    [X, base_valid] = build_eog_design( ...
        group_code, eog_ptp, eog_derivative);
    valid = base_valid & isfinite(y);
    [beta, covariance, degrees_of_freedom] = ...
        fit_ols(X(valid, :), double(y(valid)));
    contrast_labels = ["80-0"; "160-0"; "緑なし-0"; ...
        "80-緑なし"; "160-緑なし"; "80-160"; ...
        "EOG PTP +1SD"; "EOG derivative +1SD"];
    contrast_matrix = [ ...
        0  1  0  0  0  0; ...
        0  0  1  0  0  0; ...
        0  0  0  1  0  0; ...
        0  1  0 -1  0  0; ...
        0  0  1 -1  0  0; ...
        0  1 -1  0  0  0; ...
        0  0  0  0  1  0; ...
        0  0  0  0  0  1];
    n_rows = numel(contrast_labels);
    estimate = nan(n_rows, 1);
    standard_error = nan(n_rows, 1);
    if ~isempty(beta)
        estimate = contrast_matrix * beta;
        contrast_variance = diag( ...
            contrast_matrix * covariance * contrast_matrix.');
        standard_error = sqrt(max(contrast_variance, 0));
    end
    ci_lower = estimate - config.confidence_multiplier * standard_error;
    ci_upper = estimate + config.confidence_multiplier * standard_error;
    output_table = table(repmat(string(channel_name), n_rows, 1), ...
        repmat(string(band_name), n_rows, 1), ...
        repmat(string(measure_name), n_rows, 1), ...
        repmat(window_ms(1), n_rows, 1), ...
        repmat(window_ms(2), n_rows, 1), contrast_labels, ...
        estimate, standard_error, ci_lower, ci_upper, ...
        repmat(sum(valid), n_rows, 1), ...
        repmat(degrees_of_freedom, n_rows, 1), ...
        'VariableNames', {'Channel', 'Band', 'Method', ...
        'WindowStartMs', 'WindowEndMs', 'Contrast', ...
        'EstimateDb', 'StandardErrorDb', 'CI95LowerDb', ...
        'CI95UpperDb', 'ValidTrials', 'DegreesOfFreedom'});
end

function output_table = build_adjusted_timecourse_table( ...
        means, lowers, uppers, group_labels, channel_name, ...
        band_name, measure_name, times)
    output_table = table();
    for group_index = 1:numel(group_labels)
        n_times = numel(times);
        row = table(repmat(string(channel_name), n_times, 1), ...
            repmat(string(band_name), n_times, 1), ...
            repmat(string(measure_name), n_times, 1), ...
            repmat(group_labels(group_index), n_times, 1), ...
            times(:), means(group_index, :).', ...
            lowers(group_index, :).', uppers(group_index, :).', ...
            'VariableNames', {'Channel', 'Band', 'Method', ...
            'Condition', 'TimeMs', 'EOGAdjustedMeanDb', ...
            'CI95LowerDb', 'CI95UpperDb'});
        output_table = append_table_row(output_table, row);
    end
end

function create_eog_adjusted_overlay_figure( ...
        adjusted_mean, group_labels, channel_name, output_png, ...
        config, is_green_difference)
    band_indices = config.specificity_band_indices;
    measure_names = ["Total", "Induced"];
    colors = lines(4);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [30 20 1550 280 * numel(band_indices)]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, numel(band_indices), 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for band_position = 1:numel(band_indices)
        band_index = band_indices(band_position);
        for measure_index = 1:numel(measure_names)
            ax = nexttile(layout);
            if is_green_difference
                indices = 1:3;
                labels = group_labels(indices) + " - 緑なし";
                values = reshape(adjusted_mean(measure_index, ...
                    band_index, indices, :), numel(indices), []);
                green_values = reshape(adjusted_mean(measure_index, ...
                    band_index, 4, :), 1, []);
                values = values - green_values;
            else
                indices = 1:4;
                labels = group_labels;
                values = reshape(adjusted_mean(measure_index, ...
                    band_index, indices, :), numel(indices), []);
            end
            if isvector(values)
                values = reshape(values, 1, []);
            end
            for curve_index = 1:size(values, 1)
                plot_mean_line(ax, config.common_times_ms, ...
                    values(curve_index, :), colors(curve_index, :), ...
                    labels(curve_index));
            end
            xline(ax, 0, '--k', 'HandleVisibility', 'off');
            yline(ax, 0, ':k', 'HandleVisibility', 'off');
            xlim(ax, config.display_window_ms);
            grid(ax, 'on');
            xlabel(ax, '時間 (ms)');
            ylabel(ax, 'EOG調整ERSP (dB)');
            title(ax, sprintf('%s / %s', ...
                config.band_labels(band_index), ...
                measure_names(measure_index)), ...
                'Interpreter', 'none', 'Color', 'k');
            legend(ax, 'Location', 'best', 'Interpreter', 'none');
            set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
                'Box', 'on');
        end
    end
    if is_green_difference
        quantity = '条件 - 緑なし';
    else
        quantity = '全条件';
    end
    title(layout, sprintf( ...
        '%s: %s（試行EOG PTP・微分最大値で共変量調整）', ...
        channel_name, quantity), 'Interpreter', 'none', ...
        'FontWeight', 'bold', 'Color', 'k');
    exportgraphics(fig, output_png, 'Resolution', 180);
    clear cleanup;
end

function create_specificity_contrast_figure( ...
        contrast_table, channel_name, output_png, config)
    band_indices = config.specificity_band_indices;
    measure_names = ["Total", "Induced"];
    plotted_contrasts = ["80-0", "160-0", "緑なし-0", ...
        "80-緑なし", "160-緑なし", "80-160"];
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [30 20 1600 285 * numel(band_indices)]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, numel(band_indices), 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for band_position = 1:numel(band_indices)
        band_index = band_indices(band_position);
        for measure_index = 1:numel(measure_names)
            ax = nexttile(layout);
            mask = contrast_table.Band == config.band_names(band_index) & ...
                contrast_table.Method == measure_names(measure_index) & ...
                ismember(contrast_table.Contrast, plotted_contrasts);
            subset = contrast_table(mask, :);
            [~, order] = ismember(plotted_contrasts, subset.Contrast);
            valid_order = order > 0;
            subset = subset(order(valid_order), :);
            labels = plotted_contrasts(valid_order);
            x = 1:height(subset);
            lower_error = subset.EstimateDb - subset.CI95LowerDb;
            upper_error = subset.CI95UpperDb - subset.EstimateDb;
            errorbar(ax, x, subset.EstimateDb, lower_error, ...
                upper_error, 'o', 'LineWidth', 1.3, ...
                'Color', [0.10 0.35 0.70], ...
                'MarkerFaceColor', [0.10 0.35 0.70]);
            yline(ax, 0, '--k');
            xlim(ax, [0.4, max(1.6, height(subset) + 0.6)]);
            xticks(ax, x);
            xticklabels(ax, labels);
            xtickangle(ax, 25);
            grid(ax, 'on');
            ylabel(ax, '調整済み差 (dB)');
            title(ax, sprintf('%s / %s', ...
                config.band_labels(band_index), ...
                measure_names(measure_index)), ...
                'Interpreter', 'none', 'Color', 'k');
            set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
                'Box', 'on');
        end
    end
    title(layout, sprintf( ...
        '%s: EOG調整済み事前時間窓の条件差（95%% OLS CI）', ...
        channel_name), 'Interpreter', 'none', ...
        'FontWeight', 'bold', 'Color', 'k');
    exportgraphics(fig, output_png, 'Resolution', 180);
    clear cleanup;
end
