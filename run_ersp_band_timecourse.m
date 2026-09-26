function master_manifest = run_ersp_band_timecourse( ...
        input_dir, display_window_ms)
%RUN_ERSP_BAND_TIMECOURSE CWT後MATを帯域別ERSP時系列へ変換する。
%バンド帯域(没)
%
% 画面からフォルダを選ぶ:
%   run_ersp_band_timecourse
%
% フォルダとERSP平均時間範囲を指定する:
%   run_ersp_band_timecourse("C:\...\04_CWT", [-200 300])
%
% 入力フォルダ直下のF3/F4/Fz/O1/O2/Oz/PO7/PO8などを走査し、
% ファイル名や実験条件に関係なく、ersp/freqs/timesを含む全MATを処理する。
% 出力は入力フォルダと同階層の「06_バンド帯」へ保存する。

    if nargin < 1 || strlength(string(input_dir)) == 0
        selected_dir = uigetdir('', ...
            'CWT後MATのチャンネルフォルダが含まれるフォルダを選択してください');
        if isequal(selected_dir, 0)
            error('run_ersp_band_timecourse:Cancelled', ...
                '処理を中断しました。');
        end
        input_dir = string(selected_dir);
    else
        input_dir = string(input_dir);
    end

    if nargin < 2 || isempty(display_window_ms)
        display_window_ms = [-200 300];
    end
    validateattributes(display_window_ms, {'numeric'}, ...
        {'vector', 'numel', 2, 'finite', 'increasing'});
    display_window_ms = double(display_window_ms(:)');

    if ~isfolder(input_dir)
        error('run_ersp_band_timecourse:MissingInputFolder', ...
            '入力フォルダが見つかりません: %s', input_dir);
    end

    [parent_path, ~] = fileparts(input_dir);
    output_root = fullfile(parent_path, "06_バンド帯");
    if ~isfolder(output_root)
        mkdir(output_root);
    end

    band_names = [ ...
        "Delta"
        "Theta"
        "Alpha"
        "Beta"
        "LowGamma"
        "HighGamma_61_100"
        "HighGamma_101_150"
        "HighGamma_151_200"
    ];
    band_labels = [ ...
        "デルタ"
        "シータ"
        "アルファ"
        "ベータ"
        "ローガンマ"
        "ハイガンマ1"
        "ハイガンマ2"
        "ハイガンマ3"
    ];
    band_ranges = [ ...
          0.5,   3
          4,     7
          8,    13
         14,    30
         31,    60
         61,   100
        101,   150
        151,   200
    ];

    folder_list = dir(input_dir);
    folder_list = folder_list([folder_list.isdir]);
    folder_list = folder_list(~ismember({folder_list.name}, {'.', '..'}));

    master_manifest = table();
    all_band_summary = table();
    all_band_timecourses = table();
    processed_files = 0;
    skipped_files = 0;

    fprintf('\n=== ERSP帯域別時系列解析を開始 ===\n');
    fprintf('入力: %s\n', input_dir);
    fprintf('出力: %s\n', output_root);
    fprintf('表示時間範囲: %.3g–%.3g ms\n', ...
        display_window_ms(1), display_window_ms(2));

    for folder_index = 1:numel(folder_list)
        channel_name = string(folder_list(folder_index).name);
        channel_dir = string(fullfile(folder_list(folder_index).folder, ...
            folder_list(folder_index).name));
        mat_list = dir(fullfile(channel_dir, '*.mat'));

        if isempty(mat_list)
            fprintf('[スキップ] %s: MATファイルなし\n', channel_name);
            continue;
        end

        channel_output = fullfile(output_root, channel_name);
        if ~isfolder(channel_output)
            mkdir(channel_output);
        end

        for file_index = 1:numel(mat_list)
            source_file = string(fullfile(mat_list(file_index).folder, ...
                mat_list(file_index).name));
            [~, source_base] = fileparts(source_file);

            available_names = string({whos('-file', source_file).name});
            required_names = ["ersp", "freqs", "times"];
            if ~all(ismember(required_names, available_names))
                fprintf('[スキップ] %s/%s: ersp/freqs/timesなし\n', ...
                    channel_name, mat_list(file_index).name);
                skipped_files = skipped_files + 1;
                continue;
            end

            data = load(source_file, 'ersp', 'freqs', 'times');
            [ersp, freqs, times] = validate_and_orient_ersp( ...
                data, source_file);

            time_mask = times >= display_window_ms(1) & ...
                times <= display_window_ms(2);
            if ~any(time_mask)
                error('run_ersp_band_timecourse:TimeWindowUnavailable', ...
                    ['%s に指定時間範囲 %.3g–%.3g msがありません。' ...
                     '保存時間範囲: %.3g–%.3g ms'], ...
                    source_file, display_window_ms(1), ...
                    display_window_ms(2), times(1), times(end));
            end

            display_times_ms = times(time_mask);
            [band_timecourses, band_summary, timecourse_table] = ...
                calculate_band_timecourses(freqs, ...
                    ersp(:, time_mask), display_times_ms, ...
                    band_names, band_labels, band_ranges);

            safe_base = sanitize_filename(source_base);
            output_mat = fullfile(channel_output, ...
                safe_base + "_band_timecourse.mat");
            output_png = fullfile(channel_output, ...
                safe_base + "_band_timecourse.png");
            output_timecourse_csv = fullfile(channel_output, ...
                safe_base + "_band_timecourse.csv");
            output_summary_csv = fullfile(channel_output, ...
                safe_base + "_band_summary.csv");

            source_mat = source_file;
            save(output_mat, 'source_mat', 'display_window_ms', ...
                'display_times_ms', 'band_timecourses', ...
                'band_names', 'band_labels', 'band_ranges', ...
                'band_summary', '-v7');
            writetable(timecourse_table, output_timecourse_csv, ...
                'Encoding', 'UTF-8');
            writetable(band_summary, output_summary_csv, ...
                'Encoding', 'UTF-8');
            create_band_figure(display_times_ms, band_timecourses, ...
                band_summary, band_labels, band_ranges, source_base, ...
                channel_name, display_window_ms, output_png, freqs);

            current_manifest = table( ...
                channel_name, source_file, string(output_mat), ...
                string(output_png), string(output_timecourse_csv), ...
                string(output_summary_csv), ...
                min(freqs), max(freqs), numel(freqs), ...
                display_window_ms(1), display_window_ms(2), ...
                sum(band_summary.Coverage == "なし"), ...
                sum(band_summary.Coverage == "一部"), ...
                'VariableNames', { ...
                    'Channel', 'SourceMAT', 'OutputMAT', 'OutputPNG', ...
                    'OutputTimecourseCSV', 'OutputSummaryCSV', ...
                    'MinimumFrequencyHz', ...
                    'MaximumFrequencyHz', 'FrequencyBinCount', ...
                    'WindowStartMs', 'WindowEndMs', ...
                    'MissingBandCount', 'PartialBandCount'});

            if isempty(master_manifest)
                master_manifest = current_manifest;
            else
                master_manifest = [master_manifest; ...
                    current_manifest]; %#ok<AGROW>
            end

            summary_with_source = addvars(band_summary, ...
                repmat(channel_name, height(band_summary), 1), ...
                repmat(source_file, height(band_summary), 1), ...
                'Before', 1, ...
                'NewVariableNames', {'Channel', 'SourceMAT'});
            if isempty(all_band_summary)
                all_band_summary = summary_with_source;
            else
                all_band_summary = [all_band_summary; ...
                    summary_with_source]; %#ok<AGROW>
            end

            timecourse_with_source = addvars(timecourse_table, ...
                repmat(channel_name, height(timecourse_table), 1), ...
                repmat(source_file, height(timecourse_table), 1), ...
                'Before', 1, ...
                'NewVariableNames', {'Channel', 'SourceMAT'});
            if isempty(all_band_timecourses)
                all_band_timecourses = timecourse_with_source;
            else
                all_band_timecourses = [all_band_timecourses; ...
                    timecourse_with_source]; %#ok<AGROW>
            end

            processed_files = processed_files + 1;
            fprintf('[完了] %s/%s\n', ...
                channel_name, mat_list(file_index).name);
        end
    end

    if processed_files == 0
        error('run_ersp_band_timecourse:NoCwtMatFound', ...
            'ersp/freqs/timesを含むCWT後MATが見つかりませんでした。');
    end

    writetable(master_manifest, ...
        fullfile(output_root, 'ERSPバンド帯_manifest.csv'), ...
        'Encoding', 'UTF-8');
    writetable(all_band_summary, ...
        fullfile(output_root, 'ERSPバンド帯_全ファイル集計.csv'), ...
        'Encoding', 'UTF-8');
    writetable(all_band_timecourses, ...
        fullfile(output_root, 'ERSPバンド帯_全ファイル時系列.csv'), ...
        'Encoding', 'UTF-8');

    fprintf('\n=== 全工程終了 ===\n');
    fprintf('処理MAT数: %d\n', processed_files);
    fprintf('非CWT MATスキップ数: %d\n', skipped_files);
    fprintf('出力先: %s\n', output_root);

    if any(master_manifest.MissingBandCount > 0 | ...
            master_manifest.PartialBandCount > 0)
        fprintf(['注意: 入力MATに含まれない、または一部しか含まれない' ...
            '帯域があります。図とCSVのCoverage列を確認してください。\n']);
    end
end

function [ersp, freqs, times] = validate_and_orient_ersp(data, source_file)
    ersp = double(data.ersp);
    freqs = double(data.freqs(:));
    times = double(data.times(:)');

    if ~isnumeric(data.ersp) || ~isnumeric(data.freqs) || ...
            ~isnumeric(data.times)
        error('run_ersp_band_timecourse:NonNumericData', ...
            '%s のersp/freqs/timesは数値配列である必要があります。', ...
            source_file);
    end
    if size(ersp, 1) == numel(freqs) && ...
            size(ersp, 2) == numel(times)
        % 期待する向き: 周波数 x 時間
    elseif size(ersp, 2) == numel(freqs) && ...
            size(ersp, 1) == numel(times)
        ersp = ersp.';
    else
        error('run_ersp_band_timecourse:DimensionMismatch', ...
            '%s のERSPサイズとfreqs/timesが一致しません。', source_file);
    end
    if any(~isfinite(freqs)) || any(~isfinite(times))
        error('run_ersp_band_timecourse:InvalidAxes', ...
            '%s の周波数軸または時間軸にNaN/Infがあります。', source_file);
    end
    if any(diff(freqs) <= 0) || any(diff(times) <= 0)
        error('run_ersp_band_timecourse:UnsortedAxes', ...
            '%s の周波数軸または時間軸が昇順ではありません。', source_file);
    end
end

function [band_timecourses, summary, timecourse_table] = ...
        calculate_band_timecourses(freqs, ersp, display_times_ms, ...
            band_names, band_labels, band_ranges)

    n_bands = numel(band_names);
    n_times = numel(display_times_ms);
    band_timecourses = nan(n_bands, n_times);
    coverage = strings(n_bands, 1);
    frequency_bin_count = zeros(n_bands, 1);
    actual_minimum_hz = nan(n_bands, 1);
    actual_maximum_hz = nan(n_bands, 1);
    mean_ersp_db = nan(n_bands, 1);

    for band_index = 1:n_bands
        low_hz = band_ranges(band_index, 1);
        high_hz = band_ranges(band_index, 2);
        band_mask = freqs >= low_hz & freqs <= high_hz;
        frequency_bin_count(band_index) = sum(band_mask);

        if ~any(band_mask)
            coverage(band_index) = "なし";
            continue;
        end

        actual_minimum_hz(band_index) = min(freqs(band_mask));
        actual_maximum_hz(band_index) = max(freqs(band_mask));
        band_timecourses(band_index, :) = mean( ...
            ersp(band_mask, :), 1, 'omitnan');
        mean_ersp_db(band_index) = mean( ...
            band_timecourses(band_index, :), 'omitnan');

        if min(freqs) <= low_hz && max(freqs) >= high_hz
            coverage(band_index) = "全域";
        else
            coverage(band_index) = "一部";
        end
    end

    summary = table(band_names, band_labels, ...
        band_ranges(:, 1), band_ranges(:, 2), coverage, ...
        frequency_bin_count, actual_minimum_hz, actual_maximum_hz, ...
        mean_ersp_db, ...
        'VariableNames', {'Band', 'BandLabelJapanese', ...
            'DefinedMinimumHz', 'DefinedMaximumHz', 'Coverage', ...
            'FrequencyBinCount', 'ActualMinimumHz', ...
            'ActualMaximumHz', 'MeanERSPdB'});

    time_ms = repmat(display_times_ms(:), n_bands, 1);
    band = repelem(band_names, n_times);
    band_label_japanese = repelem(band_labels, n_times);
    band_coverage = repelem(coverage, n_times);
    mean_ersp_timecourse_db = reshape(band_timecourses.', [], 1);
    timecourse_table = table(time_ms, band, band_label_japanese, ...
        band_coverage, mean_ersp_timecourse_db, ...
        'VariableNames', {'TimeMs', 'Band', 'BandLabelJapanese', ...
            'Coverage', 'MeanERSPdB'});
end

function create_band_figure(display_times_ms, band_timecourses, ...
        summary, band_labels, band_ranges, source_base, channel_name, ...
        display_window_ms, output_png, freqs)

    fig = figure('Color', 'w', 'Visible', 'off', ...
        'Units', 'pixels', 'Position', [80, 50, 1350, 1200]);
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig, 4, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    finite_ersp = band_timecourses(isfinite(band_timecourses));
    if isempty(finite_ersp)
        common_ylim = [-1 1];
    else
        common_ylim = [min(finite_ersp), max(finite_ersp)];
        if common_ylim(1) == common_ylim(2)
            common_ylim = common_ylim + [-1 1];
        else
            margin = 0.08 * diff(common_ylim);
            common_ylim = common_ylim + [-margin margin];
        end
    end

    for band_index = 1:height(summary)
        ax = nexttile(layout);
        low_hz = band_ranges(band_index, 1);
        high_hz = band_ranges(band_index, 2);

        if summary.Coverage(band_index) ~= "なし"
            plot(ax, display_times_ms, ...
                band_timecourses(band_index, :), ...
                '-', 'LineWidth', 1.8, ...
                'Color', [0.05 0.35 0.70]);
            yline(ax, 0, '--k', 'LineWidth', 0.8);
            ylim(ax, common_ylim);
            grid(ax, 'on');
        else
            ylim(ax, [-1 1]);
            text(ax, mean(display_window_ms), 0, ...
                sprintf('データなし\nMAT範囲: %.3g–%.3g Hz', ...
                    min(freqs), max(freqs)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', [0.75 0.15 0.10], 'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end

        if display_window_ms(1) <= 0 && display_window_ms(2) >= 0
            xline(ax, 0, '--', ...
                'Color', [0.25 0.25 0.25], 'LineWidth', 1);
        end
        xlim(ax, display_window_ms);
        xlabel(ax, '時間 (ms)');
        ylabel(ax, '平均ERSP (dB)');
        title(ax, sprintf('%s (%.1f–%.1f Hz) [%s]', ...
            band_labels(band_index), low_hz, high_hz, ...
            summary.Coverage(band_index)), ...
            'Interpreter', 'none');
        set(ax, 'FontSize', 10, 'Box', 'on');
    end

    title(layout, sprintf( ...
        ['%s / %s: 帯域平均ERSP時系列 ' ...
         '(%.3g–%.3g ms, 0 ms = サッケード開始)'], ...
        source_base, channel_name, display_window_ms(1), ...
        display_window_ms(2)), ...
        'Interpreter', 'none', 'FontWeight', 'bold');
    exportgraphics(fig, output_png, 'Resolution', 200);
    clear cleanup;
end

function safe_name = sanitize_filename(name)
    safe_name = regexprep(string(name), '[<>:"/|?*]', '_');
    safe_name = strrep(safe_name, '\', '_');
end
