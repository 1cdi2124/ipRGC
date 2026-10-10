function master_manifest = joint_20260926( ...
        input_root, analysis_channels, ica_seed, ersp_max_db)
%RUN_JOINT_ICA_NEWTIMEF_BAND_ANALYSIS 共通ICA+従来newtimef方式のERSP解析。
%
% 実行例:
%   run_joint_ica_newtimef_band_analysis
%   run_joint_ica_newtimef_band_analysis("C:\...\03_手動サッケード")
%   run_joint_ica_newtimef_band_analysis( ...
%       "C:\...\03_手動サッケード", "Oz")
%   run_joint_ica_newtimef_band_analysis( ...
%       "C:\...\03_手動サッケード", "Oz", 10)
%   run_joint_ica_newtimef_band_analysis( ...
%       "C:\...\03_手動サッケード", [], 10)
%   run_joint_ica_newtimef_band_analysis( ...
%       "C:\...\03_手動サッケード", "Oz", 10, 1.5)
%
% 各03_手動サッケード内のSETを同一被験者の条件として扱う。
% 0Oは緑LED常灯、緑なしは無灯の眼球運動対照であり、区別する。
% 全条件へ同一フィルタ・平均基準を適用し、連結データでICAを1回だけ
% 学習する。明確なEye/Muscle/Line Noise/Channel Noiseだけを自動除去し、
% Brain確率80%未満やEOG相関だけが高いICは要確認として残す。
% 同じIC選択を全条件へ適用する。0.5-200 Hzの平均ERSPは単一の
% newtimef解析により、パワー領域でベースライン処理・試行平均した後に
% dB変換する。低周波を別手法・探索的パネルへ分離しない。
% 04ではnewtimefが返す時間・周波数グリッドを補間せず使用する。
% 05と06の条件間比較は同一軸を必須とし、不一致なら差分計算前に停止する。
%
% 出力:
%   04_ERSP_共通ICA_newtimef_全帯域           条件別ERSPとICA-QC
%   05_ERSP緑なし差分_共通ICA_newtimef_全帯域 各条件-緑なしのERSP差分
%   06_バンド帯_共通ICA_newtimef_全帯域       条件別・差分の帯域ERSP

    config = default_config();
    if nargin >= 2 && ~isempty(analysis_channels)
        analysis_channels = string(analysis_channels(:)');
        if ~all(ismember(lower(analysis_channels), ...
                lower(config.target_channels)))
            error('run_joint_ica_cwt_band_analysis:InvalidOutputChannel', ...
                '解析チャンネルは次から選んでください: %s', ...
                strjoin(config.target_channels, ', '));
        end
        config.analysis_channels = analysis_channels;
    end
    if nargin >= 3 && ~isempty(ica_seed)
        validateattributes(ica_seed, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'integer', 'nonnegative', ...
             '<=', double(intmax('uint32'))}, mfilename, 'ica_seed', 3);
        config.ica_seed = double(ica_seed);
    end
    if nargin >= 4 && ~isempty(ersp_max_db)
        validateattributes(ersp_max_db, {'numeric'}, ...
            {'scalar', 'real', 'finite', 'positive'}, ...
            mfilename, 'ersp_max_db', 4);
        config.ersp_color_limit_db = double(ersp_max_db);
    end
    if nargin < 1 || strlength(string(input_root)) == 0
        selected = uigetdir('', ...
            ['実験ルート、被験者フォルダ、または' ...
             '03_手動サッケードを選択してください']);
        if isequal(selected, 0)
            error('run_joint_ica_cwt_band_analysis:Cancelled', ...
                '処理を中断しました。');
        end
        input_root = string(selected);
    else
        input_root = string(input_root);
    end
    if ~isfolder(input_root)
        error('run_joint_ica_cwt_band_analysis:MissingInputFolder', ...
            '入力フォルダが見つかりません: %s', input_root);
    end

    initialize_eeglab(config);
    manual_folders = find_manual_saccade_folders(input_root);
    if isempty(manual_folders)
        error('run_joint_ica_cwt_band_analysis:NoInputFolder', ...
            '03_手動サッケードが見つかりません: %s', input_root);
    end

    master_manifest = table();
    fprintf('\n=== 全条件共通ICA・newtimef ERSP帯域解析を開始 ===\n');
    fprintf('検索ルート: %s\n', input_root);
    fprintf('対象フォルダ数: %d\n', numel(manual_folders));
    fprintf('解析フィルタ: %.1f–%.0f Hz / 平均基準（EOG除外）\n', ...
        config.analysis_highpass_hz, config.analysis_lowpass_hz);
    fprintf('ICA乱数シード: %u\n', config.ica_seed);
    fprintf('ERSP表示範囲: ±%.2f dB\n', config.ersp_color_limit_db);

    for folder_index = 1:numel(manual_folders)
        manual_folder = manual_folders(folder_index);
        subject_folder = string(fileparts(manual_folder));
        output_roots = make_output_roots(subject_folder, config);
        ensure_folder(output_roots.CWT);
        ensure_folder(output_roots.Difference);
        ensure_folder(output_roots.Band);
        fprintf('\n--- [%d/%d] %s ---\n', folder_index, ...
            numel(manual_folders), manual_folder);
        fprintf('  04 ERSP: %s\n', output_roots.CWT);
        fprintf('  05 緑なし差分: %s\n', output_roots.Difference);
        fprintf('  06 バンド帯: %s\n', output_roots.Band);

        try
            datasets = load_and_preprocess_conditions( ...
                manual_folder, config);
            [datasets, joint_ica_report] = run_joint_ica_cleaning( ...
                datasets, output_roots.CWT, config);
            [folder_manifest, analysis_results] = ...
                analyze_all_conditions(datasets, output_roots, ...
                    joint_ica_report, config);
            validate_all_result_axes(analysis_results, config);
            difference_manifest = create_cwt_difference_outputs( ...
                analysis_results, output_roots.Difference, config);
            if ~isempty(difference_manifest)
                writetable(difference_manifest, fullfile( ...
                    output_roots.Difference, ...
                    '共通ICA_newtimef_ERSP緑なし差分_manifest.csv'), ...
                    'Encoding', 'UTF-8');
            end
            create_all_channel_overlays(analysis_results, ...
                output_roots.Band, config);
            folder_status = "完了";
            folder_message = "";
        catch ME
            warning('run_joint_ica_cwt_band_analysis:FolderFailed', ...
                '%s の処理に失敗しました: %s', ...
                manual_folder, ME.message);
            folder_manifest = table();
            folder_status = "失敗";
            folder_message = string(getReport(ME, 'extended', ...
                'hyperlinks', 'off'));
        end

        folder_row = table(manual_folder, output_roots.CWT, ...
            output_roots.Difference, output_roots.Band, ...
            folder_status, folder_message, ...
            'VariableNames', {'ManualSaccadeFolder', 'ERSPOutputRoot', ...
                'DifferenceOutputRoot', 'BandOutputRoot', ...
                'Status', 'Message'});
        if isempty(master_manifest)
            master_manifest = folder_row;
        else
            master_manifest = [master_manifest; folder_row]; %#ok<AGROW>
        end
        if ~isempty(folder_manifest)
            writetable(folder_manifest, ...
                fullfile(output_roots.Band, '共通ICA_ERSP_manifest.csv'), ...
                'Encoding', 'UTF-8');
        end
    end

    fprintf('\n=== 全条件共通ICA・newtimef ERSP帯域解析が終了しました ===\n');
    fprintf('完了フォルダ: %d / 失敗フォルダ: %d\n', ...
        sum(master_manifest.Status == "完了"), ...
        sum(master_manifest.Status == "失敗"));
end

function config = default_config()
    config.cwt_output_relative_path = ...
        "04_ERSP_共通ICA_newtimef_全帯域";
    config.difference_output_relative_path = ...
        "05_ERSP緑なし差分_共通ICA_newtimef_全帯域";
    config.band_output_relative_path = ...
        "06_バンド帯_共通ICA_newtimef_全帯域";
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
    config.brain_confident_threshold = 0.80;
    config.artifact_remove_class_names = [ ...
        "Eye", "Muscle", "Line Noise", "Channel Noise"];
    config.artifact_remove_threshold = 0.80;
    config.eog_correlation_threshold = 0.40;

    config.display_window_ms = [-200 300];
    config.baseline_window_ms = [-250 -100];
    % 0 HzはDCで時間周波数解析できないため、最低周波数は0.5 Hzとする。
    config.newtimef_frequency_limits_hz = [0.5 200];
    config.newtimef_frequency_count = 200;
    config.newtimef_frequency_scale = 'linear';
    config.wavelet_cycles = 3;
    % 0.5 Hz・3周期（6秒）でもベースラインと表示範囲が残る長さ。
    config.epoch_window_s = [-4 4];
    config.newtimef_timesout = 400;
    config.padratio = 1;
    % 条件間比較のため固定する。空配列 [] にするとデータから自動決定。
    % ERSPヒートマップの表示範囲。解析値・保存値はクリップしない。
%-----------------------------------------------------------------------
    config.ersp_color_limit_db = 1.5;
%-----------------------------------------------------------------------
    config.mean_ersp_method = ...
        "0.5-200 Hz: epoch voltage baseline then one EEGLAB newtimef analysis (power baseline -> trial mean -> dB); no interpolation; cross-condition axes must match";
    config.sem_method = ...
        "single-trial dB values relative to the common newtimef baseline; descriptive SEM";

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
    config.band_interpretation = repmat("可視化", ...
        numel(config.band_names), 1);
    config.confidence_multiplier = 1.96;
    config.band_ersp_axis_step_db = 0.2;

    config.eeglab_fallback = ...
        "C:\Users\tatsuya\Downloads\Apps\eeglab2025.1.0";
end

function output_roots = make_output_roots(subject_folder, config)
    output_roots = struct();
    output_roots.CWT = string(fullfile(subject_folder, ...
        config.cwt_output_relative_path));
    output_roots.Difference = string(fullfile(subject_folder, ...
        config.difference_output_relative_path));
    output_roots.Band = string(fullfile(subject_folder, ...
        config.band_output_relative_path));
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
        EEG = require_finite_data(EEG);

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
    fprintf('  共通ICA: %d EEG ch, rank=%d, 学習%.0f Hz, seed=%u\n', ...
        numel(target_indices), rank_value, training.srate, config.ica_seed);
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
    report.ICASeed = config.ica_seed;
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

function EEG = require_finite_data(EEG)
    invalid = any(~isfinite(double(EEG.data)), 1);
    invalid = reshape(invalid, 1, []);
    if ~any(invalid)
        return;
    end
    error('run_joint_ica_newtimef_band_analysis:NonfiniteData', ...
        ['非有限値を含む%dサンプルがあります。線形補間は行わない設定のため、' ...
         '入力データを確認してください。'], sum(invalid));
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
    remove_mask = false(n_components, 1);
    removal_reasons = strings(n_components, 1);
    review_reasons = strings(n_components, 1);
    if isempty(classification) || ...
            size(classification, 1) ~= n_components
        error('run_joint_ica_cwt_band_analysis:InvalidICLabelOutput', ...
            'ICLabelの成分数が共通ICAの成分数と一致しません。');
    end
    for class_index = 1:numel(config.artifact_remove_class_names)
        class_name = config.artifact_remove_class_names(class_index);
        probability_column = find(strcmpi(classes, class_name), 1);
        if isempty(probability_column)
            continue;
        end
        class_probability = classification(:, probability_column);
        class_hit = isfinite(class_probability) & ...
            class_probability >= config.artifact_remove_threshold;
        remove_mask = remove_mask | class_hit;
        for component_index = find(class_hit(:))'
            removal_reasons(component_index) = append_reason( ...
                removal_reasons(component_index), sprintf( ...
                '%s=%.3f>=%.2f', class_name, ...
                class_probability(component_index), ...
                config.artifact_remove_threshold));
        end
    end

    brain_column = find(strcmpi(classes, "Brain"), 1);
    if isempty(brain_column)
        error('run_joint_ica_cwt_band_analysis:BrainClassMissing', ...
            'ICLabel出力にBrainクラスがありません。');
    end

    brain_probability = classification(:, brain_column);
    brain_uncertain = ~isfinite(brain_probability) | ...
        brain_probability < config.brain_confident_threshold;
    eog_evaluated = isfinite(eog_correlations);
    eog_hit = eog_evaluated & abs(eog_correlations) >= ...
        config.eog_correlation_threshold;
    review_mask = ~remove_mask & ...
        (brain_uncertain | eog_hit | ~eog_evaluated);

    for component_index = find(review_mask(:) & brain_uncertain(:))'
        if isfinite(brain_probability(component_index))
            reason = sprintf('Brain=%.3f<%.2f（残す）', ...
                brain_probability(component_index), ...
                config.brain_confident_threshold);
        else
            reason = 'Brain未評価（残す）';
        end
        review_reasons(component_index) = append_reason( ...
            review_reasons(component_index), reason);
    end
    for component_index = find(review_mask(:) & eog_hit(:))'
        review_reasons(component_index) = append_reason( ...
            review_reasons(component_index), sprintf( ...
            'EOG|r|=%.3f>=%.2f（単独では除去しない）', ...
            abs(eog_correlations(component_index)), ...
            config.eog_correlation_threshold));
    end
    for component_index = find(remove_mask(:) & eog_hit(:))'
        removal_reasons(component_index) = append_reason( ...
            removal_reasons(component_index), sprintf( ...
            'EOG|r|=%.3f', ...
            abs(eog_correlations(component_index))));
    end
    for component_index = find(review_mask(:) & ~eog_evaluated(:))'
        review_reasons(component_index) = append_reason( ...
            review_reasons(component_index), ...
            'EOG相関未評価（残す）');
    end
    if all(remove_mask)
        error('run_joint_ica_cwt_band_analysis:NoComponentKept', ...
            ['全ICが明確なアーチファクトとして除去対象になりました。' ...
             'ICLabel判定とICA分解を確認してください。']);
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
    ica_seed = repmat(report.ICASeed, n_components, 1);
    report_table = table(component, dominant_class, ...
        dominant_probability, report.EOGCorrelations(:), kept, review, ...
        removed, ica_seed, report.ReviewReason(:), report.RemovalReason(:), ...
        'VariableNames', {'IC', 'DominantClass', ...
        'DominantProbability', 'EOGCorrelation', 'Kept', 'Review', ...
        'Removed', 'ICASeed', 'ReviewReason', 'RemovalReason'});
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
    ica_seed = report.ICASeed;
    analysis_config = config;
    save(fullfile(qc_folder, '共通ICA_分解と除去判定.mat'), ...
        'icaweights', 'icasphere', 'icawinv', 'icachansind', ...
        'channel_labels', 'ica_seed', 'report', 'analysis_config', '-v7');

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
            decision_mark = " [採用]";
            if ismember(component_index, report.RemovedComponents)
                decision_mark = " [除去]";
            elseif ismember(component_index, report.ReviewComponents)
                decision_mark = " [要確認・採用]";
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
            ['全条件共通ICA: 採用IC = %s（要確認 = %s）/ 除去IC = %s' ...
             '（seed=%u、明確なアーチファクト >= %.0f%%だけを自動除去）'], ...
            mat2str(report.KeptComponents), ...
            mat2str(report.ReviewComponents), ...
            mat2str(report.RemovedComponents), ...
            report.ICASeed, ...
            100 * config.artifact_remove_threshold), ...
            'Interpreter', 'none');
        overall_title.Color = 'k';
        print_figure_png(fig, ...
            fullfile(qc_folder, '共通ICA_ICLabel_EOG判定.png'), 180);
        clear cleanup;
    catch ME
        warning('run_joint_ica_cwt_band_analysis:QCFigureFailed', ...
            '共通ICA QC画像を保存できませんでした: %s', ME.message);
    end
end

function [manifest, results] = analyze_all_conditions( ...
        datasets, output_roots, joint_report, config)
    n_conditions = numel(datasets);
    n_channels = numel(config.analysis_channels);
    results = repmat(struct(), n_conditions, 1);
    manifest = table();

    for condition_index = 1:n_conditions
        dataset = datasets(condition_index);
        EEG = dataset.EEG;
        EEG_epoch = epoch_for_newtimef(EEG, config);
        fprintf('  ERSP: %s（newtimef epoch %d）\n', ...
            dataset.ConditionLabel, EEG_epoch.trials);

        results(condition_index).ConditionLabel = ...
            dataset.ConditionLabel;
        results(condition_index).ConditionValue = ...
            dataset.ConditionValue;
        results(condition_index).IsGreen = dataset.IsGreen;
        results(condition_index).SourceSET = dataset.SourceSET;
        results(condition_index).BaseName = dataset.BaseName;
        channel_template = struct('Name', "", ...
            'CWTMean', [], 'CWTSEM', [], 'CWTN', [], ...
            'CWTFrequencies', [], 'CWTTimes', [], ...
            'Mean', [], ...
            'SEM', [], 'CI95Lower', [], 'CI95Upper', [], ...
            'N', [], 'Trials', []);
        results(condition_index).Channels = repmat(channel_template, ...
            n_channels, 1);

        for channel_position = 1:n_channels
            channel_name = config.analysis_channels(channel_position);
            target_position = find(strcmpi(config.target_channels, ...
                channel_name), 1);
            channel_index = dataset.TargetIndices(target_position);
            fprintf('    newtimef/帯域: %s ...\n', channel_name);
            [trial_ersp, freqs, ersp_mean, times, itc, powbase] = ...
                compute_newtimef_trial_ersp(EEG_epoch, ...
                    channel_index, config);
            [~, ersp_sem, ersp_n] = mean_sem_over_trials(trial_ersp);

            band_stats = calculate_trial_band_statistics( ...
                trial_ersp, freqs, times, config);
            band_stats = align_band_means_with_newtimef( ...
                band_stats, ersp_mean, freqs, times, config);

            cwt_channel_folder = fullfile(output_roots.CWT, channel_name);
            band_channel_folder = fullfile(output_roots.Band, channel_name);
            ensure_folder(cwt_channel_folder);
            ensure_folder(band_channel_folder);
            safe_base = safe_filename(dataset.BaseName);
            file_stem = safe_base + "_" + channel_name;
            cwt_mat = fullfile(cwt_channel_folder, ...
                file_stem + "_共通ICA_newtimef_ERSP.mat");
            cwt_png = fullfile(cwt_channel_folder, ...
                file_stem + "_共通ICA_newtimef_ERSP.png");
            band_mat = fullfile(band_channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.mat");
            band_png = fullfile(band_channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.png");
            band_csv = fullfile(band_channel_folder, ...
                file_stem + "_共通ICA_帯域ERSP.csv");

            source_set = dataset.SourceSET;
            condition_label = dataset.ConditionLabel;
            condition_value_hz = dataset.ConditionValue;
            is_green_noled = dataset.IsGreen;
            kept_components = joint_report.KeptComponents;
            review_components = joint_report.ReviewComponents;
            removed_components = joint_report.RemovedComponents;
            analysis_config = config;
            trial_count = size(trial_ersp, 3);
            save(cwt_mat, 'source_set', 'condition_label', ...
                'condition_value_hz', 'is_green_noled', ...
                'kept_components', 'review_components', ...
                'removed_components', ...
                'analysis_config', ...
                'ersp_mean', 'ersp_sem', 'ersp_n', 'freqs', ...
                'times', 'itc', 'powbase', 'trial_count', '-v7');

            band_mean = band_stats.Mean;
            band_sem = band_stats.SEM;
            band_ci95_lower = band_stats.CI95Lower;
            band_ci95_upper = band_stats.CI95Upper;
            band_n = band_stats.N;
            band_trials = band_stats.Trials;
            band_summary = band_stats.Summary;
            save(band_mat, 'source_set', 'condition_label', ...
                'condition_value_hz', 'is_green_noled', ...
                'kept_components', 'review_components', ...
                'removed_components', ...
                'analysis_config', 'times', ...
                'band_mean', 'band_sem', 'band_ci95_lower', ...
                'band_ci95_upper', 'band_n', 'band_trials', ...
                'band_summary', '-v7');
            band_table = build_band_timecourse_table( ...
                band_stats, times, dataset, channel_name, config);
            % 06のERSP CSVは、列番号ではなく変数名で不要列を削除する。
            % これにより、列タイトルと各列の値の対応を維持する。
            columns_to_remove = intersect( ...
                {'SourceSET', 'BandLabelJapanese', 'Interpretation'}, ...
                band_table.Properties.VariableNames, 'stable');
            if ~isempty(columns_to_remove)
                band_table(:, columns_to_remove) = [];
            end
            % Band列にはDelta/Theta/Alpha/...が入る。
            writetable(band_table, band_csv, 'Encoding', 'UTF-8');

            create_newtimef_figure(times, freqs, ersp_mean, ...
                dataset.ConditionLabel, channel_name, cwt_png, config, ...
                config.ersp_color_limit_db);
            create_single_condition_band_figure(times, band_stats, ...
                dataset.ConditionLabel, channel_name, band_png, config);

            channel_result = struct();
            channel_result.Name = channel_name;
            channel_result.CWTMean = ersp_mean;
            channel_result.CWTSEM = ersp_sem;
            channel_result.CWTN = ersp_n;
            channel_result.CWTFrequencies = freqs;
            channel_result.CWTTimes = times;
            channel_result.Mean = band_stats.Mean;
            channel_result.SEM = band_stats.SEM;
            channel_result.CI95Lower = band_stats.CI95Lower;
            channel_result.CI95Upper = band_stats.CI95Upper;
            channel_result.N = band_stats.N;
            channel_result.Trials = band_stats.Trials;
            results(condition_index).Channels(channel_position) = ...
                channel_result;

            row = table(dataset.SourceSET, dataset.ConditionLabel, ...
                dataset.ConditionValue, dataset.IsGreen, channel_name, ...
                dataset.RawSaccadeEvents, dataset.DuplicatesRemoved, ...
                EEG_epoch.trials, config.ica_seed, ...
                config.ersp_color_limit_db, ...
                string(mat2str(joint_report.KeptComponents)), ...
                string(mat2str(joint_report.ReviewComponents)), ...
                string(mat2str(joint_report.RemovedComponents)), ...
                string(cwt_mat), string(cwt_png), string(band_mat), ...
                string(band_png), string(band_csv), ...
                'VariableNames', {'SourceSET', 'Condition', ...
                'ConditionValueHz', 'IsGreenNoLED', 'Channel', ...
                'RawSaccadeEvents', 'DuplicateEventsRemoved', ...
                'NewtimefTrials', 'ICASeed', 'ERSPDisplayMaxDb', ...
                'CommonKeptICs', 'CommonReviewICs', ...
                'CommonRemovedICs', 'ERSPMAT', 'ERSPPNG', ...
                'BandMAT', 'BandPNG', 'BandCSV'});
            manifest = append_table_row(manifest, row);

            clear trial_ersp band_trials
        end
    end
end

function EEG_epoch = epoch_for_newtimef(EEG, config)
    EEG_epoch = pop_epoch(EEG, {char(config.event_type)}, ...
        config.epoch_window_s, 'epochinfo', 'yes');
    if EEG_epoch.trials < 2
        error('run_joint_ica_newtimef_band_analysis:TooFewEpochs', ...
            '端部を除外するとnewtimef用epochが2試行未満です。');
    end
    % 従来のICA_CWT_ERS_0702.mと同じく、時間領域のベースラインを
    % 除去してからnewtimef側でもパワーベースライン補正を行う。
    EEG_epoch = pop_rmbase(EEG_epoch, config.baseline_window_ms);
    EEG_epoch = eeg_checkset(EEG_epoch, 'eventconsistency');
end

function [trial_ersp, freqs, mean_ersp, times, itc, powbase] = ...
        compute_newtimef_trial_ersp(EEG_epoch, channel_index, config)
    signal = double(reshape(EEG_epoch.data(channel_index, :, :), ...
        EEG_epoch.pnts, EEG_epoch.trials));
    [mean_ersp, itc, powbase, times, freqs, ~, ~, tf] = ...
        newtimef(signal, EEG_epoch.pnts, ...
        [EEG_epoch.xmin EEG_epoch.xmax] * 1000, EEG_epoch.srate, ...
        config.wavelet_cycles, 'plotersp', 'off', ...
        'plotitc', 'off', 'plotphase', 'off', ...
        'baseline', config.baseline_window_ms, ...
        'freqs', config.newtimef_frequency_limits_hz, ...
        'nfreqs', config.newtimef_frequency_count, ...
        'freqscale', config.newtimef_frequency_scale, ...
        'timesout', config.newtimef_timesout, ...
        'padratio', config.padratio, 'verbose', 'off');

    freqs = double(freqs(:));
    times = double(times(:)');
    mean_ersp = double(real(mean_ersp));
    itc = double(itc);
    powbase = double(powbase);
    tf = normalize_tf_dimensions(tf, numel(freqs), ...
        EEG_epoch.trials);
    linear_power = max(abs(tf).^2, realmin('double'));
    baseline_power = 10 .^ (powbase(:) / 10);
    if numel(baseline_power) ~= numel(freqs) || ...
            any(~isfinite(baseline_power) | baseline_power <= 0)
        error('run_joint_ica_newtimef_band_analysis:InvalidBaseline', ...
            'newtimefが返したベースラインが不正です。');
    end
    % trialbase='off' のnewtimefと同じ、全試行共通の線形パワー
    % ベースラインを用いる。平均ERSPそのものはnewtimef出力を採用する。
    trial_ersp = 10 * log10(linear_power ./ ...
        reshape(baseline_power, [], 1, 1));
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
        trial_ersp, freqs, times, config)
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
        band_mask = freqs >= low_hz & freqs <= high_hz;
        excluded_line_bin_count(band_index) = 0;
        frequency_bin_count(band_index) = sum(band_mask);
        if ~any(band_mask)
            coverage(band_index) = "なし";
            stats.Trials{band_index} = nan(n_times, 0);
            continue;
        end

        actual_minimum_hz(band_index) = min(freqs(band_mask));
        actual_maximum_hz(band_index) = max(freqs(band_mask));
        if min(freqs) <= low_hz && max(freqs) >= high_hz
            coverage(band_index) = "全域";
        else
            coverage(band_index) = "一部";
        end
        values = mean(trial_ersp(band_mask, :, :), 1, 'omitnan');
        values = reshape(values, n_times, size(trial_ersp, 3));
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

function stats = align_band_means_with_newtimef( ...
        stats, ersp_mean, freqs, times, config)
    if size(ersp_mean, 2) ~= numel(times)
        error('run_joint_ica_newtimef_band_analysis:BandTimeMismatch', ...
            '帯域平均用ERSPの時間軸が一致しません。');
    end
    for band_index = 1:numel(config.band_names)
        band_mask = freqs >= config.band_ranges(band_index, 1) & ...
            freqs <= config.band_ranges(band_index, 2);
        if ~any(band_mask)
            continue;
        end
        stats.Mean(band_index, :) = mean( ...
            ersp_mean(band_mask, :), 1, 'omitnan');
        stats.CI95Lower(band_index, :) = ...
            stats.Mean(band_index, :) - ...
            config.confidence_multiplier * stats.SEM(band_index, :);
        stats.CI95Upper(band_index, :) = ...
            stats.Mean(band_index, :) + ...
            config.confidence_multiplier * stats.SEM(band_index, :);
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

function create_newtimef_figure(times, freqs, ersp_mean, ...
        condition_label, channel_name, output_png, config, color_limit)
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [80 60 1250 760]);
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);
    imagesc(ax, times, freqs, ersp_mean);
    set(ax, 'YDir', 'normal', 'FontSize', 11, 'Box', 'on', ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    colormap(ax, turbo);
    xlim(ax, config.display_window_ms);
    ylim(ax, [min(freqs), max(freqs)]);
    hold(ax, 'on');
    xline(ax, 0, '--m', 'LineWidth', 1.4);
    xlabel(ax, '時間 (ms)');
    ylabel(ax, '周波数 (Hz)');

    if isempty(color_limit)
        finite_values = ersp_mean(isfinite(ersp_mean));
        color_limit = finite_percentile(abs(finite_values), 98);
    end
    if isfinite(color_limit) && color_limit > 0
        clim(ax, [-color_limit color_limit]);
    end
    colorbar_handle = colorbar(ax);
    colorbar_handle.Label.String = '試行平均ERSP (dB)';
    colorbar_handle.Color = 'k';
    title(ax, sprintf( ...
        '%s / %s: 共通ICA後newtimef ERSP（0.5–200 Hz・無補間）', ...
        condition_label, channel_name), ...
        'Interpreter', 'none', 'FontWeight', 'bold');
    print_figure_png(fig, output_png, 200);
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

    % 条件内の全帯域で共通の縦軸を使用する。まず全帯域をまとめて
    % 表示余白込みの範囲を求め、その外側を0.2 dB単位に丸める。
    common_y_limits = curve_limits( ...
        stats.Mean, zeros(size(stats.Mean)));
    axis_step_db = config.band_ersp_axis_step_db;
    common_y_limits(1) = axis_step_db * ...
        floor(common_y_limits(1) / axis_step_db);
    common_y_limits(2) = axis_step_db * ...
        ceil(common_y_limits(2) / axis_step_db);
    if diff(common_y_limits) < axis_step_db
        common_y_limits = common_y_limits + ...
            [-axis_step_db, axis_step_db];
    end

    for band_index = 1:numel(config.band_names)
        ax = nexttile(layout);
        mean_values = stats.Mean(band_index, :);
        plot_mean_line(ax, times, mean_values, line_color, "");
        hold(ax, 'on');
        xline(ax, 0, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1);
        yline(ax, 0, ':k', 'LineWidth', 0.8);
        xlim(ax, config.display_window_ms);
        ylim(ax, common_y_limits);
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
    print_figure_png(fig, output_png, 200);
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

function print_figure_png(fig, output_png, resolution)
    % exportgraphicsを多数回呼び出した際のgraphics handshaking timeoutを
    % 避けるため、非表示figureは従来型のprintでPNGへ保存する。
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    print(fig, char(output_png), '-dpng', ...
        sprintf('-r%d', resolution));
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

function manifest = create_cwt_difference_outputs( ...
        results, output_root, config)
    manifest = table();
    if isempty(results)
        return;
    end

    is_green = reshape([results.IsGreen], [], 1);
    green_indices = find(is_green);
    if isempty(green_indices)
        fprintf('  緑なし条件がないためERSP差分を省略します。\n');
        return;
    end
    if numel(green_indices) > 1
        warning('run_joint_ica_cwt_band_analysis:MultipleGreenCWT', ...
            '緑なし条件が複数あるためERSP差分を一意に作成できません。');
        return;
    end

    green_index = green_indices(1);
    condition_indices = find(~is_green);
    for channel_position = 1:numel(config.analysis_channels)
        channel_name = config.analysis_channels(channel_position);
        channel_folder = fullfile(output_root, channel_name);
        ensure_folder(channel_folder);
        green_channel = results(green_index).Channels(channel_position);

        for difference_index = 1:numel(condition_indices)
            condition_index = condition_indices(difference_index);
            condition_channel = ...
                results(condition_index).Channels(channel_position);
            % 05と06の両方で事前検証した同一軸上で差分を計算する。
            difference_ersp_mean = condition_channel.CWTMean - ...
                green_channel.CWTMean;
            difference_ersp_sem = sqrt( ...
                condition_channel.CWTSEM .^ 2 + ...
                green_channel.CWTSEM .^ 2);
            condition_n = condition_channel.CWTN;
            green_n = green_channel.CWTN;
            freqs = green_channel.CWTFrequencies;
            times = green_channel.CWTTimes;

            % 旧run_ersp_band_timecourseが読み込める互換変数名。
            ersp = difference_ersp_mean;
            condition_source_set = results(condition_index).SourceSET;
            green_source_set = results(green_index).SourceSET;
            condition_label = results(condition_index).ConditionLabel;
            difference_label = condition_label + " - 緑なし";
            analysis_config = config;
            metadata = struct( ...
                'method', ...
                "Common-ICA contrast; one 0.5-200 Hz conventional newtimef mean", ...
                'equation', "condition ERSP mean - green-no-LED ERSP mean", ...
                'sem_method', "sqrt(condition SEM^2 + green SEM^2)", ...
                'interpolation', "none; matching condition axes required", ...
                'interpretation_note', ...
                "newtimef ERSP後の条件差であり、生波形から眼球運動を物理的に除去したものではない");

            safe_base = safe_filename(results(condition_index).BaseName);
            file_stem = safe_base + "_" + channel_name + ...
                "_共通ICA_newtimef_ERSP_緑なし差分";
            output_mat = fullfile(channel_folder, file_stem + ".mat");
            output_png = fullfile(channel_folder, file_stem + ".png");
            save(output_mat, 'ersp', 'difference_ersp_mean', ...
                'difference_ersp_sem', 'condition_n', 'green_n', ...
                'freqs', 'times', 'condition_source_set', ...
                'green_source_set', 'condition_label', ...
                'difference_label', 'metadata', 'analysis_config', '-v7');
            create_newtimef_figure(times, freqs, difference_ersp_mean, ...
                difference_label, channel_name, output_png, config, ...
                config.ersp_color_limit_db);

            row = table(string(channel_name), condition_label, ...
                condition_source_set, green_source_set, ...
                string(output_mat), string(output_png), ...
                min(freqs), max(freqs), times(1), times(end), ...
                'VariableNames', {'Channel', 'DifferenceCondition', ...
                    'ConditionSourceSET', 'GreenNoLEDSourceSET', ...
                    'OutputMAT', 'ERSPPNG', 'MinimumFrequencyHz', ...
                    'MaximumFrequencyHz', 'TimeStartMs', 'TimeEndMs'});
            manifest = append_table_row(manifest, row);
        end
    end
end

function validate_all_result_axes(results, config)
%VALIDATE_ALL_RESULT_AXES 05と06の全条件を同一グリッドに限定する。
    if isempty(results)
        return;
    end
    reference_index = find(reshape([results.IsGreen], [], 1), 1);
    if isempty(reference_index)
        reference_index = 1;
    end
    for channel_position = 1:numel(config.analysis_channels)
        channel_name = config.analysis_channels(channel_position);
        reference_channel = ...
            results(reference_index).Channels(channel_position);
        for condition_index = 1:numel(results)
            condition_channel = ...
                results(condition_index).Channels(channel_position);
            validate_cwt_difference_axes(condition_channel, ...
                reference_channel, ...
                results(condition_index).ConditionLabel, ...
                results(reference_index).ConditionLabel, channel_name);
        end
    end
end

function validate_cwt_difference_axes(condition_channel, reference_channel, ...
        condition_label, reference_label, channel_name)
    same_frequencies = isequal(size(condition_channel.CWTFrequencies), ...
        size(reference_channel.CWTFrequencies)) && all(abs( ...
        condition_channel.CWTFrequencies(:) - ...
        reference_channel.CWTFrequencies(:)) < 1e-9);
    same_times = isequal(size(condition_channel.CWTTimes), ...
        size(reference_channel.CWTTimes)) && all(abs( ...
        condition_channel.CWTTimes(:) - ...
        reference_channel.CWTTimes(:)) < 1e-9);
    same_data_size = isequal(size(condition_channel.CWTMean), ...
        size(reference_channel.CWTMean));
    if ~(same_frequencies && same_times && same_data_size)
        error('run_joint_ica_newtimef_band_analysis:ERSPAxisMismatch', ...
            '%s / %sのERSP軸が基準条件%sと一致しません。', ...
            condition_label, channel_name, reference_label);
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
        times = results(1).Channels(channel_position).CWTTimes;
        [raw_mean, raw_sem, raw_n] = collect_channel_arrays( ...
            results, channel_position, config);
        raw_png = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.png');
        raw_mat = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.mat');
        raw_csv = fullfile(channel_folder, ...
            '全条件_帯域ERSP重ね描き_SEM.csv');
        create_overlay_figure(times, raw_mean, ...
            condition_labels, channel_name, raw_png, ...
            config, false);
        raw_table = build_overlay_table(times, ...
            raw_mean, raw_sem, raw_n, condition_labels, ...
            source_sets, channel_name, config);
        writetable(raw_table, raw_csv, 'Encoding', 'UTF-8');
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
            numel(config.band_names), numel(times));
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
        create_overlay_figure(times, ...
            difference_mean, difference_labels, ...
            channel_name, difference_png, config, true);
        difference_table = build_difference_overlay_table( ...
            times, difference_mean, difference_sem, ...
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
    n_times = size(results(1).Channels(channel_position).Mean, 2);
    array_size = [n_conditions, numel(config.band_names), ...
        n_times];
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
    print_figure_png(fig, output_png, 200);
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
