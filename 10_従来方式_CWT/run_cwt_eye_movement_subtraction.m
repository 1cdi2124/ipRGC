function master_manifest = run_cwt_eye_movement_subtraction(input_dir)
%RUN_CWT_EYE_MOVEMENT_SUBTRACTION CWTフォルダ内の全チャンネルを一括補正する。
%緑なし差分
%
% 画面から選択する場合:
%   run_cwt_eye_movement_subtraction
%
% フォルダを直接指定する場合:
%   run_cwt_eye_movement_subtraction("C:\...\04_CWT")
%
% 選択フォルダ直下のチャンネルフォルダ（F3、Ozなど）を走査し、
% ファイル名に「緑なし」を含むMATに対応する0O/80O/160Oを自動検出する。

    if nargin < 1 || strlength(string(input_dir)) == 0
        selected_dir = uigetdir('', ...
            'CWT結果のチャンネルフォルダが含まれるフォルダを選択してください');
        if isequal(selected_dir, 0)
            error('run_cwt_eye_movement_subtraction:Cancelled', ...
                '処理を中断しました。');
        end
        input_dir = string(selected_dir);
    else
        input_dir = string(input_dir);
    end

    if ~isfolder(input_dir)
        error('run_cwt_eye_movement_subtraction:MissingInputFolder', ...
            '入力フォルダが見つかりません: %s', input_dir);
    end

    % ICA_CWT_ERS_0702.mと同様に、入力フォルダと同じ階層へ出力する。
    [parent_path, ~] = fileparts(input_dir);
    output_root = fullfile(parent_path, "05_CWT緑なし差分");
    if ~isfolder(output_root)
        mkdir(output_root);
    end

    condition_tokens = ["0O", "80O", "160O"];
    condition_names = ["0", "80", "160"];

    folder_list = dir(input_dir);
    folder_list = folder_list([folder_list.isdir]);
    folder_list = folder_list(~ismember({folder_list.name}, {'.', '..'}));

    master_manifest = table();
    processed_groups = 0;

    fprintf('\n=== CWT眼球運動寄与補正を開始 ===\n');
    fprintf('入力: %s\n', input_dir);
    fprintf('出力: %s\n', output_root);

    for folder_index = 1:numel(folder_list)
        channel_name = string(folder_list(folder_index).name);
        channel_dir = string(fullfile(folder_list(folder_index).folder, ...
            folder_list(folder_index).name));

        mat_list = dir(fullfile(channel_dir, '*.mat'));
        if isempty(mat_list)
            fprintf('[スキップ] %s: MATファイルなし\n', channel_name);
            continue;
        end

        mat_names = string({mat_list.name});
        mat_keys = strings(size(mat_names));
        for name_index = 1:numel(mat_names)
            mat_keys(name_index) = make_pairing_key(mat_names(name_index));
        end

        noled_indices = find(contains(mat_names, "緑なし"));
        if isempty(noled_indices)
            fprintf('[スキップ] %s: 「緑なし」を含むMATなし\n', channel_name);
            continue;
        end
        if numel(unique(mat_keys(noled_indices))) ~= numel(noled_indices)
            error('run_cwt_eye_movement_subtraction:DuplicateNoLed', ...
                ['%s に同じ測定を表す「緑なし」MATが複数あります。\n' ...
                 '対象: %s'], channel_name, ...
                strjoin(mat_names(noled_indices), ', '));
        end

        channel_output = fullfile(output_root, channel_name);
        for noled_index = noled_indices
            noled_name = mat_names(noled_index);
            noled_file = fullfile(channel_dir, noled_name);
            noled_key = mat_keys(noled_index);

            led_files = strings(numel(condition_tokens), 1);
            for condition_index = 1:numel(condition_tokens)
                token_mask = has_condition_token( ...
                    mat_names, condition_tokens(condition_index));
                candidate_indices = find(token_mask & mat_keys == noled_key);

                if isempty(candidate_indices)
                    error('run_cwt_eye_movement_subtraction:MissingPair', ...
                        ['%s の対応ファイルがありません。\n' ...
                         '緑なし: %s\n正規化キー: %s\n不足条件: %s'], ...
                        channel_name, noled_name, noled_key, ...
                        condition_tokens(condition_index));
                end
                if numel(candidate_indices) > 1
                    error('run_cwt_eye_movement_subtraction:AmbiguousPair', ...
                        ['%s の条件%sに候補が複数あります。\n' ...
                         '緑なし: %s\n候補: %s'], ...
                        channel_name, condition_tokens(condition_index), ...
                        noled_name, ...
                        strjoin(mat_names(candidate_indices), ', '));
                end
                led_files(condition_index) = fullfile( ...
                    channel_dir, mat_names(candidate_indices));
            end

            fprintf('\n--- %s: %s を基準に処理 ---\n', ...
                channel_name, noled_name);
            current_manifest = subtract_cwt_eye_movement( ...
                led_files, noled_file, channel_output, condition_names);

            group_name = erase(noled_name, ".mat");
            current_manifest.channel = repmat(channel_name, ...
                height(current_manifest), 1);
            current_manifest.group = repmat(group_name, ...
                height(current_manifest), 1);
            current_manifest = movevars(current_manifest, ...
                {'channel', 'group'}, 'Before', 1);

            if isempty(master_manifest)
                master_manifest = current_manifest;
            else
                master_manifest = [master_manifest; current_manifest]; %#ok<AGROW>
            end
            processed_groups = processed_groups + 1;
        end
    end

    if processed_groups == 0
        error('run_cwt_eye_movement_subtraction:NoPairsFound', ...
            ['処理可能な組が見つかりませんでした。入力フォルダ直下に、' ...
             'チャンネル別フォルダと「*_緑なし_*.mat」があるか確認してください。']);
    end

    master_manifest_path = fullfile(output_root, ...
        'CWT眼球運動成分補正_全体manifest.csv');
    writetable(master_manifest, master_manifest_path, 'Encoding', 'UTF-8');

    fprintf('\n=== 全工程終了 ===\n');
    fprintf('処理グループ数: %d\n', processed_groups);
    fprintf('出力先: %s\n', output_root);
    fprintf('全体manifest: %s\n', master_manifest_path);
end

function key = make_pairing_key(file_name)
% 条件名、測定日、Windowsのコピー番号を除去して対応キーを作る。
% 例:
%   kojima_緑なし_0618_F3.mat -> kojima_f3
%   kojima_80O_0617_F3.mat    -> kojima_f3
%   takahira_緑なし_0527 (2)_Oz.mat -> takahira_oz

    [~, base_name] = fileparts(file_name);
    base_name = regexprep(base_name, '緑なし', '');
    base_name = regexprep(base_name, ...
        '(^|_)(0O|80O|160O)(_|$)', '$1$3');
    base_name = regexprep(base_name, '\s*\(\d+\)', '');
    base_name = regexprep(base_name, '\s*（\d+）', '');
    base_name = regexprep(base_name, ...
        '(^|_)\d{4}(_|$)', '$1$2');
    base_name = regexprep(base_name, '[_\s-]+', '_');
    base_name = regexprep(base_name, '^_+|_+$', '');
    key = lower(string(base_name));
end

function mask = has_condition_token(file_names, condition_token)
% 0Oが80Oや160Oの一部として誤検出されないよう境界付きで検索する。

    mask = false(size(file_names));
    token = regexptranslate('escape', char(condition_token));
    pattern = ['(^|_)' token '(_|$)'];
    for name_index = 1:numel(file_names)
        [~, base_name] = fileparts(file_names(name_index));
        mask(name_index) = ~isempty( ...
            regexp(char(base_name), pattern, 'once'));
    end
end

function manifest = subtract_cwt_eye_movement( ...
        led_files, noled_file, output_dir, condition_names)
% 1組のLED条件と緑なし条件を検証し、CWT差分を保存するローカル関数。

    led_files = string(led_files);
    noled_file = string(noled_file);
    output_dir = string(output_dir);
    condition_names = string(condition_names);

    if numel(condition_names) ~= numel(led_files)
        error('run_cwt_eye_movement_subtraction:ConditionCountMismatch', ...
            '条件名とLED条件ファイルの数が一致しません。');
    end
    if ~isfile(noled_file)
        error('run_cwt_eye_movement_subtraction:MissingNoLedFile', ...
            '緑なし条件が見つかりません: %s', noled_file);
    end
    if any(~isfile(led_files))
        error('run_cwt_eye_movement_subtraction:MissingLedFile', ...
            'LED条件ファイルが見つかりません。');
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    required_fields = {'ersp', 'itc', 'ers', 'powbase', 'times', 'freqs'};
    noled = load(noled_file);
    validate_cwt_data(noled, required_fields, noled_file);

    n_conditions = numel(led_files);
    corrected = repmat(struct( ...
        'ersp', [], 'itc', [], 'itc_magnitude', [], 'ers', [], ...
        'powbase', [], 'times', [], 'freqs', [], 'metadata', struct()), ...
        n_conditions, 1);

    all_ersp = [];
    all_itc_magnitude = [];
    all_ers = [];

    for condition_index = 1:n_conditions
        led = load(led_files(condition_index));
        validate_cwt_data(led, required_fields, led_files(condition_index));
        validate_matching_axes( ...
            led, noled, led_files(condition_index), noled_file);

        corrected(condition_index).ersp = led.ersp - noled.ersp;
        corrected(condition_index).itc = led.itc - noled.itc;
        corrected(condition_index).itc_magnitude = ...
            abs(led.itc) - abs(noled.itc);
        corrected(condition_index).ers = led.ers - noled.ers;
        corrected(condition_index).powbase = led.powbase - noled.powbase;
        corrected(condition_index).times = led.times;
        corrected(condition_index).freqs = led.freqs;
        corrected(condition_index).metadata = struct( ...
            'condition', condition_names(condition_index), ...
            'led_source', led_files(condition_index), ...
            'noled_source', noled_file, ...
            'created_at', string(datetime('now', ...
                'Format', 'yyyy-MM-dd HH:mm:ss Z')), ...
            'method', "CWT-domain condition subtraction", ...
            'ersp_equation', "LED.ersp - NOLED.ersp", ...
            'ers_equation', "LED.ers - NOLED.ers", ...
            'itc_equation', "LED.itc - NOLED.itc", ...
            'itc_magnitude_equation', ...
                "abs(LED.itc) - abs(NOLED.itc)", ...
            'interpretation_note', ...
                "CWT後の条件差分であり、生波形から眼球運動成分を物理的に削除したものではない");

        all_ersp = [all_ersp; ...
            corrected(condition_index).ersp(:)]; %#ok<AGROW>
        all_itc_magnitude = [all_itc_magnitude; ...
            corrected(condition_index).itc_magnitude(:)]; %#ok<AGROW>
        all_ers = [all_ers; ...
            corrected(condition_index).ers(:)]; %#ok<AGROW>
    end

    common_limits = [symmetric_limit(all_ersp), ...
        symmetric_limit(all_itc_magnitude), symmetric_limit(all_ers)];

    condition = strings(n_conditions, 1);
    led_source = strings(n_conditions, 1);
    noled_source = repmat(noled_file, n_conditions, 1);
    output_mat = strings(n_conditions, 1);
    output_png = strings(n_conditions, 1);
    n_freqs = zeros(n_conditions, 1);
    n_times = zeros(n_conditions, 1);
    time_start_ms = zeros(n_conditions, 1);
    time_end_ms = zeros(n_conditions, 1);
    freq_start_hz = zeros(n_conditions, 1);
    freq_end_hz = zeros(n_conditions, 1);

    for condition_index = 1:n_conditions
        condition(condition_index) = condition_names(condition_index);
        led_source(condition_index) = led_files(condition_index);

        [~, led_base_name] = fileparts(led_files(condition_index));
        safe_name = sanitize_filename(led_base_name);
        mat_path = fullfile(output_dir, ...
            safe_name + "_眼球運動補正.mat");
        png_path = fullfile(output_dir, ...
            safe_name + "_眼球運動補正.png");

        ersp = corrected(condition_index).ersp;
        itc = corrected(condition_index).itc;
        itc_magnitude = corrected(condition_index).itc_magnitude;
        ers = corrected(condition_index).ers;
        powbase = corrected(condition_index).powbase;
        times = corrected(condition_index).times;
        freqs = corrected(condition_index).freqs;
        metadata = corrected(condition_index).metadata;

        save(mat_path, 'ersp', 'itc', 'itc_magnitude', 'ers', ...
            'powbase', 'times', 'freqs', 'metadata', '-v7');
        create_difference_figure(corrected(condition_index), ...
            condition_names(condition_index), common_limits, png_path);

        output_mat(condition_index) = mat_path;
        output_png(condition_index) = png_path;
        n_freqs(condition_index) = numel(freqs);
        n_times(condition_index) = numel(times);
        time_start_ms(condition_index) = times(1);
        time_end_ms(condition_index) = times(end);
        freq_start_hz(condition_index) = freqs(1);
        freq_end_hz(condition_index) = freqs(end);

        fprintf('[完了] %s -> %s\n', ...
            condition_names(condition_index), mat_path);
    end

    manifest = table(condition, led_source, noled_source, ...
        output_mat, output_png, n_freqs, n_times, time_start_ms, ...
        time_end_ms, freq_start_hz, freq_end_hz);

    [~, noled_base_name] = fileparts(noled_file);
    manifest_name = "CWT補正_manifest_" + ...
        sanitize_filename(noled_base_name) + ".csv";
    writetable(manifest, fullfile(output_dir, manifest_name), ...
        'Encoding', 'UTF-8');
end

function validate_cwt_data(data, required_fields, file_path)
    missing = required_fields(~isfield(data, required_fields));
    if ~isempty(missing)
        error('run_cwt_eye_movement_subtraction:MissingVariables', ...
            '%s に必要な変数がありません: %s', ...
            file_path, strjoin(missing, ', '));
    end

    for field_index = 1:numel(required_fields)
        field_name = required_fields{field_index};
        value = data.(field_name);
        if ~isnumeric(value)
            error('run_cwt_eye_movement_subtraction:NonNumericVariable', ...
                '%s の %s が数値配列ではありません。', ...
                file_path, field_name);
        end
        if any(~isfinite(value(:)))
            error('run_cwt_eye_movement_subtraction:NonFiniteValue', ...
                '%s の %s にNaNまたはInfがあります。', ...
                file_path, field_name);
        end
    end
end

function validate_matching_axes(led, noled, led_file, noled_file)
    matrix_fields = {'ersp', 'itc', 'ers'};
    for field_index = 1:numel(matrix_fields)
        field_name = matrix_fields{field_index};
        if ~isequal(size(led.(field_name)), size(noled.(field_name)))
            error('run_cwt_eye_movement_subtraction:SizeMismatch', ...
                '%s の %s サイズが緑なし条件と一致しません。', ...
                led_file, field_name);
        end
    end
    if ~isequal(size(led.powbase), size(noled.powbase))
        error('run_cwt_eye_movement_subtraction:PowbaseSizeMismatch', ...
            '%s のpowbaseサイズが緑なし条件と一致しません。', led_file);
    end

    tolerance = 1e-10;
    if ~isequal(size(led.times), size(noled.times)) || ...
            any(abs(led.times(:) - noled.times(:)) > tolerance)
        error('run_cwt_eye_movement_subtraction:TimeAxisMismatch', ...
            '%s と %s の時間軸が一致しません。', led_file, noled_file);
    end
    if ~isequal(size(led.freqs), size(noled.freqs)) || ...
            any(abs(led.freqs(:) - noled.freqs(:)) > tolerance)
        error('run_cwt_eye_movement_subtraction:FrequencyAxisMismatch', ...
            '%s と %s の周波数軸が一致しません。', led_file, noled_file);
    end
end

function value = symmetric_limit(data)
    data = data(isfinite(data));
    if isempty(data)
        value = 1;
        return;
    end
    value = max(abs(data));
    if value == 0
        value = 1;
    end
end

function safe_name = sanitize_filename(name)
    safe_name = regexprep(string(name), '[<>:"/|?*]', '_');
    safe_name = strrep(safe_name, '\', '_');
end

function create_difference_figure( ...
        data, condition_name, common_limits, output_path)

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'pixels', 'Position', [100, 100, 1050, 1050]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 3, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    targets = {data.ersp, data.itc_magnitude, data.ers};
    titles = {'ERSP差分 (dB)', 'ITC振幅差分', 'ERS差分 (dB)'};
    colorbar_labels = { ...
        'LED - 緑なし (dB)', 'LED - 緑なし', 'LED - 緑なし (dB)'};

    for plot_index = 1:3
        ax = nexttile(layout);
        imagesc(ax, data.times, data.freqs, targets{plot_index});
        axis(ax, 'xy');
        xline(ax, 0, '--k', 'LineWidth', 1);
        colormap(ax, jet(256));
        clim(ax, [-common_limits(plot_index), ...
            common_limits(plot_index)]);
        cb = colorbar(ax);
        cb.Label.String = colorbar_labels{plot_index};
        ylabel(ax, '周波数 (Hz)');
        title(ax, titles{plot_index});
        set(ax, 'FontSize', 11, 'Box', 'on');
    end

    xlabel(layout, '時間 (ms)');
    title(layout, "CWT眼球運動寄与補正: " + condition_name, ...
        'Interpreter', 'none', 'FontWeight', 'bold');
    exportgraphics(fig, output_path, 'Resolution', 200);
    clear cleanup;
end
