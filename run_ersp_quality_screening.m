function [cell_qc, participant_qc, output_dir] = run_ersp_quality_screening( ...
        input_root, output_root, options)
%RUN_ERSP_QUALITY_SCREENING 既存newtimef ERSPの非破壊QC＋三要因ANOVA用CSV。
% 実行例:
%   run_ersp_quality_screening
%   run_ersp_quality_screening("C:\...\01_CWT_joint_new_NoGreen_Brain")
%   run_ersp_quality_screening("C:\...\matsumoto_0818")
%   [cells, subjects, out] = run_ersp_quality_screening("C:\...", "D:\QC")
%   run_ersp_quality_screening("C:\...", "", struct('min_trials_review', 30))
%   run_ersp_quality_screening("C:\...", "", struct('mode',"anova_only"))
%   run_ersp_quality_screening("C:\...", "", struct('mode',"qc_only"))
%
% 引数（[] / "" は既定値）:
% 1 input_root : 被験者一覧、被験者、03_手動サッケード、04_ERSPまたはその
%               06帯域またはそのチャンネルフォルダ。既定: uigetdir。直下の被験者だけを探索し、
%               異なる解析版・実験を再帰的に混ぜない。
%               04/06の表示上限別出力と従来の接尾辞なし出力を自動検出する。
%               複数版なら選択画面を表示。1版なら自動選択。同じ版の04/06だけを読む。
%               04/06またはそのchを直接選んだ場合、その版を使用する。
% 2 output_root: 出力の親フォルダ。既定: 選択した範囲の08_ERSP_QC。
%               毎回新規のrun_日時フォルダを作成し、前回出力を上書きしない。
% 3 options    : default_config()と同名の項目を持つscalar struct。
%               例: struct('display_limit_db',1.5,'min_trials_review',30)。
%               未知の項目・不正な閾値はエラー。既定値は本ファイル内に記載。
% options.mode: all（既定: QC＋07）/ qc_only / anova_only（既存04/06から07だけ）/
%               wide_only（既存ANOVA_all_cells_QC.csvから07を整形）。
% ANOVAは0～300 msの既存06帯域dB平均。06がない条件のみ保存済み04から復元する。
% 選択した表示上限をQC画像にも適用。metric_frequency_hzとANOVA帯域は変更しない。
% optionsでフォルダ名/表示範囲を明示した場合は、その指定を優先する。
% QCの候補判定では自動除外しない。
% 07も毎回run_日時を新規作成する。anova_only/wide_onlyの第3出力は07の保存先。
%
% CWT_20261006等の04出力 *_共通ICA_newtimef_ERSP.mat を使用する。
% ersp_mean, times, freqsを必須とし、trial_count, ersp_n,
% condition_label, source_set, analysis_configも読み込む。
% EEGLAB / Statistics Toolbox不要。MATLAB R2020b以降を想定。
%
% 重要:
% - 試作の形状指標・閾値。ノイズ原因や電極抵抗を確定する分類器ではない。
% - 採用候補は「本QCで警告なし」であり、品質の保証ではない。
% - 赤い面積、条件差、群の期待パターン、p値で判定しない。
% - 形状による候補判定は元解析のベースライン区間を使用する。
%   表示区間全体の指標も記録するが、イベント後の反応だけでは格上げしない。
% - 初版はSET/PNG画像の画素を読まず、試行の除去・ICA再計算を行わない。
% - MAT/SET/PNG/既存ANOVAを変更・削除しない。新規07に三要因用CSVを保存する。
% - 本来の活動も警告され得る。教授の判定例と照合して閾値を固定してから使用。
%
% 出力: QC_cells.csv / QC_participants.csv / QC_review.html / QC_results.mat
%       画像: images（採用候補も含む全セル）、boards（警告被験者の全条件・ch一覧）
%       boardsのファイル名は被験者フォルダ名.png（例: hasegawa_0807_F.png）。
%       HTMLの切替: 全て / 両方（要確認＋除外候補）/ 要確認 / 除外候補。
% 判定: AcceptCandidate=採用候補 / Review=要確認 / ExcludeCandidate=除外候補。
% 欠測・重複・計算不能を「採用候補」にしない。被験者IDはフォルダ名そのまま。

    program_file = string(mfilename('fullpath')) + ".m";
    program_source = fileread(program_file); % 判定コードも実行開始時のものを記録
    config = default_config();
    explicit = struct('ERSPFolder',false,'BandFolder',false,'DisplayFrequency',false);
    if nargin >= 3 && ~isempty(options)
        if ~isstruct(options) || ~isscalar(options)
            error('erspQC:InvalidOptions', 'optionsはscalar structで指定してください。');
        end
        names = fieldnames(options);
        explicit.ERSPFolder = isfield(options,'ersp_folder_name');
        explicit.BandFolder = isfield(options,'band_folder_name');
        explicit.DisplayFrequency = isfield(options,'display_frequency_hz');
        for k = 1:numel(names)
            if ~isfield(config, names{k})
                error('erspQC:UnknownOption', '未知の設定: %s', names{k});
            end
            config.(names{k}) = options.(names{k});
        end
    end
    config = validate_config(config);
    if nargin < 1 || isempty(input_root) || strlength(string(input_root)) == 0
        selected = uigetdir('', '被験者一覧、被験者、04/06または既存cells表のフォルダを選択');
        if isequal(selected, 0)
            error('erspQC:Cancelled', '選択を中断しました。');
        end
        input_root = string(selected);
    end
    input_root = scalar_folder(input_root);
    custom_output = nargin >= 2 && ~isempty(output_root) && strlength(string(output_root)) > 0;
    if config.mode == "wide_only"
        if ~custom_output, output_root = ""; end
        [output_dir,anova_info] = export_existing_anova(input_root,string(output_root),config);
        cell_qc=table(); participant_qc=table(); %#ok<NASGU>
        save(fullfile(output_dir,'ANOVA_export_results.mat'),'anova_info','config','program_source','program_file');
        return;
    end
    [config, analysis_variant] = resolve_analysis_variant(input_root,config,explicit);
    config = validate_config(config);
    fprintf('使用する解析版: %s\n  04: %s\n  06: %s\n', ...
        analysis_variant.Label,config.ersp_folder_name,config.band_folder_name);
    [subjects, scope_root] = discover_subjects(input_root, config);
    if nargin < 2 || isempty(output_root) || strlength(string(output_root)) == 0
        output_root = fullfile(scope_root, '08_ERSP_QC');
    end
    output_root = string(output_root);
    if ~isscalar(output_root) || ismissing(output_root) || strlength(output_root) == 0
        error('erspQC:InvalidOutput', '出力先は1つのフォルダで指定してください。');
    end
    if strlength(config.anova_output_root)==0
        if custom_output, stats_parent=output_root; else, stats_parent=scope_root; end
        config.anova_output_root=string(fullfile(stats_parent,config.anova_folder_name));
    end
    % 同名セル/被験者は黙って統合しない。hasegawaの2フォルダも別ID。
    if numel(unique(subjects.ParticipantID)) ~= height(subjects)
        error('erspQC:DuplicateParticipantID', '同名被験者があります。解析版を1つに絞って選択してください。');
    end
    if config.mode == "anova_only"
        output_dir=unique_output_dir(config.anova_output_root);
        anova_info=export_band_anova(subjects,table(),output_dir,config);
        cell_qc=table(); participant_qc=table(); %#ok<NASGU>
        save(fullfile(output_dir,'ANOVA_export_results.mat'),'anova_info','config','program_source','program_file','analysis_variant');
        return;
    end
    rows = repmat(empty_record(), 0, 1);
    fprintf('ERSP QC試作版: %d名。元データは変更しません。\n', height(subjects));
    for s = 1:height(subjects)
        fprintf('  [%d/%d] %s\n', s, height(subjects), subjects.ParticipantID(s));
        for c = 1:numel(config.channels)
            channel_dir = fullfile(subjects.ERSPFolder(s), config.channels(c));
            if strlength(subjects.ERSPFolder(s)) == 0 || ~isfolder(channel_dir)
                continue;
            end
            files = dir(fullfile(channel_dir, '*_共通ICA_newtimef_ERSP.mat'));
            [~, order] = sort(string({files.name}));
            files = files(order);
            for f = 1:numel(files)
                row = empty_record();
                row.ParticipantID = subjects.ParticipantID(s);
                row.Channel = config.channels(c);
                row.SourceMAT = string(fullfile(files(f).folder, files(f).name));
                row.SourcePNG = replace(row.SourceMAT, '.mat', '.png');
                row.FileID = string(sprintf('cell_%05d', numel(rows) + 1));
                row = inspect_file(row, config);
                % 未知条件も残してReviewにする。既知条件の意図的な絞り込みのみ除く。
                if ismember(row.Condition, ["0 Hz", "80 Hz", "160 Hz", "NoGreen"]) && ...
                        ~ismember(row.Condition, config.conditions)
                    continue;
                end
                rows(end + 1, 1) = row; %#ok<AGROW>
            end
        end
    end
    if isempty(rows)
        error('erspQC:NoERSPFiles', ...
            ['選択した解析版のERSP MATが0件です。入力範囲と04フォルダを確認してください。' ...
             '\n04: %s\n06のみを集計する場合はmode="anova_only"を使用してください。'], ...
            config.ersp_folder_name);
    else
        cell_qc = struct2table(rows);
        cell_qc = flag_duplicates_and_channel_outliers(cell_qc, config);
    end
    participant_qc = summarize_subjects(subjects, cell_qc, config);
    % スキャンと判定後に初めて出力。mkdir失敗は入力へ影響しない。
    output_dir = unique_output_dir(output_root);
    anova_info=struct('Status',"Disabled",'OutputDir',"",'Files',strings(0,1));
    if config.mode == "all"
        try
            stats_dir=unique_output_dir(config.anova_output_root);
            anova_info=export_band_anova(subjects,participant_qc,stats_dir,config);
            save(fullfile(stats_dir,'ANOVA_export_results.mat'),'anova_info','config','program_source','program_file','analysis_variant');
        catch err
            anova_info=struct('Status',"Failed",'OutputDir',"",'Files',strings(0,1), ...
                'Error',string(getReport(err,'extended','hyperlinks','off')));
            warning('erspQC:ANOVAExportFailed','07 ANOVA用CSVの出力に失敗（QCは継続）: %s',err.message);
        end
    end
    mkdir(fullfile(output_dir, 'images'));
    mkdir(fullfile(output_dir, 'boards'));
    flagged = find(cell_qc.CandidateCode ~= "AcceptCandidate");
    render_indices = (1:height(cell_qc))';
    fprintf('指標計算終了: %dセル（警告%dセル）。全セルの画像を作成します。\n', ...
        height(cell_qc),numel(flagged));
    for k = 1:numel(render_indices)
        if k == 1 || mod(k,20) == 0 || k == numel(render_indices)
            fprintf('  セル画像 [%d/%d]\n',k,numel(render_indices));
        end
        i = render_indices(k);
        image_rel = "images/" + cell_qc.FileID(i) + ".png";
        try
            render_cell(cell_qc(i, :), fullfile(output_dir, image_rel), config);
            cell_qc.ReviewImage(i) = image_rel;
        catch err
            cell_qc.RenderError(i) = string(err.message);
            warning('erspQC:RenderFailure', '画像作成失敗 %s: %s', cell_qc.FileID(i), err.message);
        end
    end
    if config.make_subject_boards
        for s = 1:height(participant_qc)
            if participant_qc.CandidateCode(s) == "AcceptCandidate"
                continue;
            end
            fprintf('  被験者一覧画像: %s\n',participant_qc.ParticipantID(s));
            board_rel = "boards/" + participant_qc.ParticipantID(s) + ".png";
            try
                render_subject_board(participant_qc.ParticipantID(s), cell_qc, ...
                    fullfile(output_dir, board_rel), config);
                participant_qc.ReviewBoard(s) = board_rel;
            catch err
                participant_qc.RenderError(s) = string(err.message);
                warning('erspQC:BoardFailure', '一覧画像作成失敗 %s: %s', ...
                    participant_qc.ParticipantID(s), err.message);
            end
        end
    end
    writetable(cell_qc, fullfile(output_dir, 'QC_cells.csv'), 'Encoding', 'UTF-8');
    writetable(participant_qc, fullfile(output_dir, 'QC_participants.csv'), 'Encoding', 'UTF-8');
    html_path = fullfile(output_dir, 'QC_review.html');
    write_html(html_path, cell_qc, participant_qc, config, anova_info);
    qc_config = config;
    save(fullfile(output_dir, 'QC_results.mat'), 'cell_qc', 'participant_qc', ...
        'qc_config', 'input_root', 'subjects', 'program_file', 'program_source', 'anova_info', 'analysis_variant', '-v7');
    fprintf('完了: %dセル / %d名 / 警告%dセル\n%s\n', ...
        height(cell_qc), height(participant_qc), numel(flagged), html_path);
    fprintf('試作閾値による候補です。QC候補での自動除外・既存ANOVAの変更はしていません。\n');
    if config.open_report
        web(char(html_path), '-browser');
    end
end

function config = default_config()
% ここで閾値を変更できる。数値は検証用の仮設定であり標準基準ではない。
    config.mode="all"; % all / qc_only / anova_only / wide_only
    config.band_folder_name="06_バンド帯_共通ICA_newtimef_全帯域";
    config.anova_folder_name="07_ANOVA_バンド帯_チャンネル";
    config.anova_output_root=""; % 空なら入力範囲（保存先指定時は指定先）配下の07
    config.anova_time_window_ms=[0 300];
    config.anova_band_names=["Delta","Theta","Alpha","Beta","LowGamma", ...
        "HighGamma_61_100","HighGamma_101_150","HighGamma_151_200"];
    config.channels = ["F3", "F4", "Fz", "O1", "O2", "Oz", "PO7", "PO8"];
    config.conditions = ["0 Hz", "80 Hz", "160 Hz", "NoGreen"];
    config.ersp_folder_name = "04_ERSP_共通ICA_newtimef_全帯域";
    config.fallback_baseline_ms = [-250 -100]; % 元MATに設定がないときのみ使用
    config.display_window_ms = [-300 300]; % baseline全体も見える共通表示
    config.display_frequency_hz = [0.5 200];
    config.display_limit_db = 3.0; % 表示だけ。指標の値をクリップしない
    config.metric_frequency_hz = [30 200]; % 低周波の時間にじみを縦筋判定へ使わない
    config.temporal_flank_ms = [10 30]; % 中心から左右10-30msの中央値と比較
    config.frequency_flank_hz = [3 10]; % 中心から上下3-10Hzの中央値と比較
    config.min_baseline_points = 5;
    config.min_trials_review = 20; % 仮値。必要試行数は研究に合わせて決める
    config.min_trials_exclude = 2; % 2未満を除外候補（元解析も2以上を要求）
    config.vertical_review_db = 0.75;
    config.vertical_severe_db = 2.0;
    config.horizontal_review_db = 0.75;
    config.horizontal_severe_db = 2.0;
    config.roughness_review_db_per_10ms = 1.0;
    config.roughness_severe_db_per_10ms = 3.0;
    config.severe_metrics_for_exclude = 2; % 複数の強いbaseline異常で除外候補
    config.channel_outlier_z = 3.5; % 同一被験者・条件のch間。単独ではReviewのみ
    config.channel_outlier_ratio = 3.0;
    config.channel_log_mad_floor = 0.05; % ch指標が同値付近でもZを暴走させない
    config.min_channels_for_outlier = 4;
    config.subject_exclude_cell_fraction = 0.25; % 25%が除外候補なら被験者も候補
    config.make_subject_boards = true;
    config.open_report = false;
end

function config = validate_config(config)
    config.mode=string(config.mode);
    if ~isscalar(config.mode) || ~ismember(config.mode,["all","qc_only","anova_only","wide_only"])
        error('erspQC:InvalidMode','modeはall/qc_only/anova_only/wide_onlyを指定してください。');
    end
    validateattributes(config.anova_time_window_ms,{'numeric'},{'real','finite','numel',2,'increasing'});
    config.anova_time_window_ms=double(config.anova_time_window_ms(:)');
    config.anova_band_names=string(config.anova_band_names(:)');
    if isempty(config.anova_band_names) || any(ismissing(config.anova_band_names)) || ...
            numel(unique(config.anova_band_names))~=numel(config.anova_band_names) || ...
            ~all(cellfun(@isvarname,cellstr(config.anova_band_names)))
        error('erspQC:InvalidBands','anova_band_namesは重複のない有効な帯域名で指定してください。');
    end
    config.anova_output_root=string(config.anova_output_root);
    if ~isscalar(config.anova_output_root) || ismissing(config.anova_output_root)
        error('erspQC:InvalidANOVAOutput','anova_output_rootは単一のパスで指定してください。');
    end
    for field={'band_folder_name','anova_folder_name'}
        name=string(config.(field{1}));
        if ~isscalar(name) || ismissing(name) || isempty(regexp(char(name),'^[^/\\]+$','once')) || ...
                ismember(name,[".",".."])
            error('erspQC:InvalidStageName','%sは単一フォルダ名で指定してください。',field{1});
        end
        config.(field{1})=name;
    end
    config.channels = string(config.channels(:)');
    config.conditions = string(config.conditions(:)');
    if isempty(config.channels) || isempty(config.conditions) || ...
            any(ismissing(config.channels)) || any(ismissing(config.conditions)) || ...
            numel(unique(config.channels)) ~= numel(config.channels) || ...
            numel(unique(config.conditions)) ~= numel(config.conditions) || ...
            ~all(ismember(config.conditions, ["0 Hz", "80 Hz", "160 Hz", "NoGreen"])) || ...
            ~all(ismember(config.channels, ["F3", "F4", "Fz", "O1", "O2", "Oz", "PO7", "PO8"]))
        error('erspQC:InvalidLevels', '条件・chは既定の候補から重複なく指定してください。');
    end
    windows = {'fallback_baseline_ms', 'display_window_ms', ...
        'display_frequency_hz', 'metric_frequency_hz', ...
        'temporal_flank_ms', 'frequency_flank_hz'};
    for k = 1:numel(windows)
        x = config.(windows{k});
        validateattributes(x, {'numeric'}, {'real', 'finite', 'vector', 'numel', 2});
        if x(2) <= x(1)
            error('erspQC:InvalidWindow', '%sは昇順の2値で指定してください。', windows{k});
        end
        config.(windows{k}) = double(x(:)');
    end
    if config.fallback_baseline_ms(2) > 0 || config.temporal_flank_ms(1) <= 0 || ...
            config.frequency_flank_hz(1) <= 0 || config.metric_frequency_hz(1) <= 0
        error('erspQC:InvalidWindow', 'baselineはイベント前、周辺窓・周波数は正値にしてください。');
    end
    positive = {'display_limit_db','vertical_review_db','vertical_severe_db', ...
        'horizontal_review_db','horizontal_severe_db', ...
        'roughness_review_db_per_10ms','roughness_severe_db_per_10ms', ...
        'channel_outlier_z','channel_outlier_ratio','channel_log_mad_floor'};
    for k = 1:numel(positive)
        validateattributes(config.(positive{k}), {'numeric'}, {'scalar','real','finite','positive'});
    end
    integer_fields = {'min_baseline_points','min_trials_review','min_trials_exclude', ...
        'severe_metrics_for_exclude','min_channels_for_outlier'};
    for k = 1:numel(integer_fields)
        validateattributes(config.(integer_fields{k}), {'numeric'}, ...
            {'scalar','real','finite','integer','positive'});
    end
    if config.severe_metrics_for_exclude > 3 || config.min_baseline_points < 3 || ...
            config.min_trials_review < config.min_trials_exclude || ...
            config.vertical_severe_db < config.vertical_review_db || ...
            config.horizontal_severe_db < config.horizontal_review_db || ...
            config.roughness_severe_db_per_10ms < config.roughness_review_db_per_10ms
        error('erspQC:InvalidThreshold', 'review/severe閾値または最小数の関係が不正です。');
    end
    validateattributes(config.subject_exclude_cell_fraction, {'numeric'}, ...
        {'scalar','real','finite','>',0,'<=',1});
    validateattributes(config.make_subject_boards, {'logical'}, {'scalar'});
    validateattributes(config.open_report, {'logical'}, {'scalar'});
    config.ersp_folder_name = string(config.ersp_folder_name);
    if ~isscalar(config.ersp_folder_name) || ismissing(config.ersp_folder_name) || ...
            isempty(regexp(char(config.ersp_folder_name), '^04_[^/\\]+$', 'once'))
        error('erspQC:InvalidStageName', 'ersp_folder_nameは04_で始まる単一フォルダ名にしてください。');
    end
end

function root = scalar_folder(value)
    root = string(value);
    if ~isscalar(root) || ismissing(root) || ~isfolder(root)
        error('erspQC:MissingFolder', '入力フォルダがありません。');
    end
    root = string(char(java.io.File(char(root)).getCanonicalPath()));
end

function [config, selected] = resolve_analysis_variant(root,config,explicit)
% 被験者直下だけを探索。表示上限の異なる04/06や従来名を混ぜない。
    ersp_base = strip_display_suffix(config.ersp_folder_name);
    band_base = strip_display_suffix(config.band_folder_name);
    [parent,name] = fileparts(root);
    [grandparent,parent_name] = fileparts(parent);
    [direct_suffix,~,direct] = parse_analysis_stage(name,ersp_base,band_base);
    if direct
        folders = string(parent);
    else
        [direct_suffix,~,direct] = parse_analysis_stage(parent_name,ersp_base,band_base);
        if direct
            folders = string(grandparent);
        elseif startsWith(string(name),"03_")
            folders = string(parent);
        elseif has_analysis_stage(root,ersp_base,band_base) || ...
                isfolder(fullfile(root,'03_手動サッケード'))
            folders = string(root);
        else
            entries = dir(root);
            entries = entries([entries.isdir] & ...
                ~startsWith(string({entries.name}),[".","07_","08_"]));
            folders = strings(0,1);
            for index = 1:numel(entries)
                folder = string(fullfile(root,entries(index).name));
                if has_analysis_stage(folder,ersp_base,band_base) || ...
                        isfolder(fullfile(folder,'03_手動サッケード'))
                    folders(end+1,1) = folder; %#ok<AGROW>
                end
            end
        end
    end
    candidates = collect_analysis_variants(folders,ersp_base,band_base);
    if explicit.ERSPFolder || explicit.BandFolder
        % 旧options形式も維持。片方だけ指定した場合は同じ接尾辞の相方を使う。
        [suffix,maximum,matched] = parse_analysis_stage(config.ersp_folder_name,ersp_base,band_base);
        if ~explicit.ERSPFolder
            [suffix,maximum,matched] = parse_analysis_stage(config.band_folder_name,ersp_base,band_base);
            if matched, config.ersp_folder_name = ersp_base+suffix; end
        elseif ~explicit.BandFolder && matched
            config.band_folder_name = band_base+suffix;
        end
        [band_suffix,~,band_matched] = parse_analysis_stage(config.band_folder_name,ersp_base,band_base);
        if matched && band_matched && suffix~=band_suffix
            error('erspQC:VariantMismatch','04と06の表示上限が異なります。同じ版を指定してください。');
        end
        method = "ExplicitOptions";
    else
        if config.mode ~= "anova_only"
            candidates = candidates(candidates.ERSPSubjectCount>0,:);
        end
        if isempty(candidates)
            error('erspQC:NoAnalysisOutputs', ...
                '対象範囲に使用可能な04/06出力がありません: %s',root);
        end
        if direct
            choice = find(candidates.Suffix==direct_suffix,1);
            if isempty(choice)
                error('erspQC:NoERSPFiles','選択した版に04出力がありません。06のみならanova_onlyを使用してください。');
            end
            method = "DirectFolder";
        elseif height(candidates)==1
            choice = 1; method = "OnlyVariant";
        else
            labels = strings(height(candidates),1);
            for index = 1:height(candidates)
                labels(index) = variant_label(candidates.Suffix(index)) + ...
                    sprintf('   [04: %d名 / 06: %d名]', ...
                    candidates.ERSPSubjectCount(index),candidates.BandSubjectCount(index));
            end
            [choice,confirmed] = listdlg('ListString',cellstr(labels), ...
                'SelectionMode','single','InitialValue',1, ...
                'Name','ERSP解析版の選択','ListSize',[700 260], ...
                'PromptString',{'使用する解析版を1つ選択してください。', ...
                    '同じ版の04と06を使用します。別の表示上限の出力は混ぜません。'}, ...
                'OKString','この版で実行','CancelString','中断');
            if ~confirmed
                error('erspQC:VariantSelectionCancelled','解析版の選択を中断しました。出力は作成していません。');
            end
            method = "Dialog";
        end
        suffix = candidates.Suffix(choice);
        maximum = candidates.DisplayMaxHz(choice);
        config.ersp_folder_name = candidates.ERSPFolderName(choice);
        config.band_folder_name = candidates.BandFolderName(choice);
    end
    if ~explicit.DisplayFrequency && isfinite(maximum)
        config.display_frequency_hz(2) = maximum;
    end
    selected = struct('Label',variant_label(suffix),'Suffix',suffix, ...
        'DisplayMaxHz',maximum,'SelectionMethod',method, ...
        'ERSPFolderName',config.ersp_folder_name,'BandFolderName',config.band_folder_name, ...
        'Candidates',candidates);
end

function base = strip_display_suffix(name)
    base = string(regexprep(char(name),'_表示[0-9]+(?:\.[0-9]+)?Hz以下$',''));
end

function [suffix,maximum,matched,kind] = parse_analysis_stage(name,ersp_base,band_base)
    suffix = ""; maximum = NaN; matched = false; kind = 0;
    bases = [ersp_base,band_base];
    for index = 1:2
        if string(name)==bases(index)
            matched = true; kind = index; return;
        end
        pattern = ['^' regexptranslate('escape',char(bases(index))) ...
            '_表示([0-9]+(?:\.[0-9]+)?)Hz以下$'];
        tokens = regexp(char(name),pattern,'tokens','once');
        if ~isempty(tokens)
            maximum = str2double(tokens{1});
            if isfinite(maximum) && maximum>0.5
                suffix = extractAfter(string(name),strlength(bases(index)));
                matched = true; kind = index; return;
            end
        end
    end
end

function found = has_analysis_stage(folder,ersp_base,band_base)
    entries = dir(folder); entries = entries([entries.isdir]);
    found = false;
    for index = 1:numel(entries)
        [~,~,matched] = parse_analysis_stage(entries(index).name,ersp_base,band_base);
        if matched, found = true; return; end
    end
end

function candidates = collect_analysis_variants(folders,ersp_base,band_base)
    candidates = table('Size',[0 6], ...
        'VariableTypes',{'string','double','string','string','double','double'}, ...
        'VariableNames',{'Suffix','DisplayMaxHz','ERSPFolderName','BandFolderName', ...
        'ERSPSubjectCount','BandSubjectCount'});
    for folder = folders(:)'
        entries = dir(folder); entries = entries([entries.isdir]);
        for index = 1:numel(entries)
            [suffix,maximum,matched,kind] = parse_analysis_stage(entries(index).name,ersp_base,band_base);
            if ~matched, continue; end
            row = find(candidates.Suffix==suffix,1);
            if isempty(row)
                candidates = [candidates; table(suffix,maximum,ersp_base+suffix,band_base+suffix,0,0, ...
                    'VariableNames',candidates.Properties.VariableNames)]; %#ok<AGROW>
                row = height(candidates);
            end
            if kind==1
                candidates.ERSPSubjectCount(row) = candidates.ERSPSubjectCount(row)+1;
            else
                candidates.BandSubjectCount(row) = candidates.BandSubjectCount(row)+1;
            end
        end
    end
    candidates = sortrows(candidates,'DisplayMaxHz');
end

function label = variant_label(suffix)
    if strlength(suffix)==0
        label = "従来名（表示上限の接尾辞なし）";
    else
        label = extractAfter(suffix,1);
    end
end

function [subjects, scope_root] = discover_subjects(root, config)
    [parent, name] = fileparts(root);
    [grandparent, parent_name] = fileparts(parent);
    if ismember(string(name),[config.ersp_folder_name,config.band_folder_name])
        scope_root = string(parent);
    elseif ismember(string(parent_name),[config.ersp_folder_name,config.band_folder_name])
        scope_root = string(grandparent); % chを選択してもその被験者の全ch
    elseif startsWith(string(name), "03_")
        scope_root = string(parent);
    else
        scope_root = root;
    end
    candidate = fullfile(scope_root, config.ersp_folder_name);
    if isfolder(candidate) || isfolder(fullfile(scope_root,config.band_folder_name))
        if ~isfolder(candidate), candidate=""; end
        [~, id] = fileparts(scope_root);
        subjects = table(string(id), string(scope_root), string(candidate), ...
            'VariableNames', {'ParticipantID','SubjectFolder','ERSPFolder'});
        return;
    end
    if isfolder(fullfile(scope_root,'03_手動サッケード'))
        [~, id] = fileparts(scope_root);
        subjects = table(string(id),string(scope_root),"", ...
            'VariableNames', {'ParticipantID','SubjectFolder','ERSPFolder'});
        return;
    end
    entries = dir(scope_root);
    entries = entries([entries.isdir] & ~ismember(string({entries.name}), [".",".."]));
    [~, order] = sort(string({entries.name}));
    entries = entries(order);
    ids = strings(0,1); folders = strings(0,1); stages = strings(0,1);
    for k = 1:numel(entries)
        if startsWith(string(entries(k).name), [".","07_","08_"])
            continue;
        end
        folder = fullfile(scope_root, entries(k).name);
        stage = fullfile(folder, config.ersp_folder_name);
        if isfolder(stage) || isfolder(fullfile(folder,config.band_folder_name))
            if ~isfolder(stage), stage=""; end
            ids(end+1,1) = string(entries(k).name); %#ok<AGROW>
            folders(end+1,1) = string(folder); %#ok<AGROW>
            stages(end+1,1) = string(stage); %#ok<AGROW>
        elseif isfolder(fullfile(folder, '03_手動サッケード'))
            % 対象被験者だが04が未作成。欠測として残す。
            ids(end+1,1) = string(entries(k).name); %#ok<AGROW>
            folders(end+1,1) = string(folder); %#ok<AGROW>
            stages(end+1,1) = ""; %#ok<AGROW>
        end
    end
    if isempty(ids)
        error('erspQC:NoSubjects', ...
            '04_ERSPがありません。被験者または被験者一覧フォルダを選択してください: %s', root);
    end
    subjects = table(ids, folders, stages, ...
        'VariableNames', {'ParticipantID','SubjectFolder','ERSPFolder'});
end

function row = empty_record()
    row = struct('FileID',"", 'ParticipantID',"", 'Condition',"Unknown", ...
        'OriginalCondition',"", 'Channel',"", 'TrialCount',NaN, ...
        'MinimumStoredN',NaN, 'BaselineStartMs',NaN, 'BaselineEndMs',NaN, ...
        'BaselinePointCount',0, 'FrequencyCount',0, 'TimePointCount',0, ...
        'MedianTimeStepMs',NaN, 'NonfiniteFraction',NaN, ...
        'BaselineVerticalDb',NaN, 'BaselineHorizontalDb',NaN, ...
        'BaselineRoughnessDbPer10Ms',NaN, 'DisplayVerticalDb',NaN, ...
        'DisplayHorizontalDb',NaN, 'DisplayClippedFraction',NaN, ...
        'ChannelRoughnessRatio',NaN, 'ChannelRobustZ',NaN, ...
        'ReviewMetricCount',0, 'SevereMetricCount',0, 'CandidateCode',"Review", ...
        'CandidateLabel',"要確認", 'Reasons',"", 'SourceMAT',"", ...
        'SourcePNG',"", 'ReviewImage',"", 'RenderError',"");
end

function row = inspect_file(row, config)
    try
        data = load_available(row.SourceMAT, {'ersp_mean','times','freqs','trial_count', ...
            'ersp_n','condition_label','source_set','analysis_config'});
        if isfield(data, 'condition_label')
            row.OriginalCondition = scalar_text(data.condition_label);
        end
        row.Condition = normalize_condition(row.OriginalCondition);
        [~, base] = fileparts(row.SourceMAT);
        if row.Condition == "Unknown"
            row.Condition = normalize_condition(string(base));
        end
        if row.Condition == "Unknown"
            row = add_reason(row, "UnknownCondition");
        end
        [map, times, freqs] = checked_map(data);
        row.TimePointCount = numel(times); row.FrequencyCount = numel(freqs);
        row.MedianTimeStepMs = median(diff(times));
        row.NonfiniteFraction = mean(~isfinite(map(:)));
        if any(~isfinite(map(:)))
            row = set_candidate(row, "ExcludeCandidate");
            row = add_reason(row, "NonfiniteERSP");
        else
            row = set_candidate(row, "AcceptCandidate");
        end
        if row.Condition == "Unknown"
            row = raise_candidate(row, "Review");
        end
        if isfield(data, 'trial_count') && isnumeric(data.trial_count) && ...
                isscalar(data.trial_count) && isfinite(data.trial_count) && ...
                data.trial_count >= 0 && fix(data.trial_count) == data.trial_count
            row.TrialCount = double(data.trial_count);
            if row.TrialCount < config.min_trials_exclude
                row = raise_candidate(row, "ExcludeCandidate");
                row = add_reason(row, "TooFewTrialsForAnalysis");
            elseif row.TrialCount < config.min_trials_review
                row = raise_candidate(row, "Review");
                row = add_reason(row, "LowTrialCount_Provisional");
            end
        else
            row = raise_candidate(row, "Review");
            row = add_reason(row, "MissingOrInvalidTrialCount");
        end
        if isfield(data, 'ersp_n')
            n = data.ersp_n;
            if isnumeric(n) && isequal(size(n),size(map)) && isreal(n) && ...
                    all(isfinite(n(:))) && all(n(:) >= 0)
                row.MinimumStoredN = min(n(:));
                if row.MinimumStoredN < config.min_trials_exclude
                    row = raise_candidate(row, "ExcludeCandidate");
                    row = add_reason(row, "TooFewStoredN");
                elseif row.MinimumStoredN < config.min_trials_review
                    row = raise_candidate(row, "Review");
                    row = add_reason(row, "LowStoredN_Provisional");
                end
            else
                row = raise_candidate(row, "Review");
                row = add_reason(row, "InvalidStoredN");
            end
        end
        baseline = config.fallback_baseline_ms;
        if isfield(data,'analysis_config') && isstruct(data.analysis_config) && ...
                isfield(data.analysis_config,'baseline_window_ms')
            baseline = double(data.analysis_config.baseline_window_ms);
        else
            row = raise_candidate(row,"Review");
            row = add_reason(row,"BaselineSettingMissing_UsingFallback");
        end
        if numel(baseline) ~= 2 || any(~isfinite(baseline)) || ...
                baseline(2) <= baseline(1) || baseline(2) > 0
            error('erspQC:InvalidBaseline', '元MATのbaseline設定が不正です。');
        end
        row.BaselineStartMs = baseline(1); row.BaselineEndMs = baseline(2);
        baseline_mask = times >= baseline(1) & times <= baseline(2);
        display_mask = times >= config.display_window_ms(1) & times <= config.display_window_ms(2);
        freq_mask = freqs >= config.metric_frequency_hz(1) & freqs <= config.metric_frequency_hz(2);
        row.BaselinePointCount = sum(baseline_mask);
        if sum(display_mask) < config.min_baseline_points
            row = raise_candidate(row,"Review");
            row = add_reason(row,"InsufficientDisplayCoverage");
        end
        if times(1) > baseline(1) || times(end) < baseline(2)
            row = raise_candidate(row,"Review");
            row = add_reason(row,"PartialBaselineCoverage");
        end
        if row.BaselinePointCount < config.min_baseline_points || sum(freq_mask) < 5
            row = raise_candidate(row,"Review");
            row = add_reason(row,"InsufficientMetricCoverage");
            return;
        end
        % 形状判定ではbaseline内のデータだけ使用し、イベント後を混ぜない。
        [row.BaselineVerticalDb, row.BaselineHorizontalDb, ...
            row.BaselineRoughnessDbPer10Ms] = shape_metrics( ...
            map(freq_mask,baseline_mask), times(baseline_mask), freqs(freq_mask), config);
        if sum(display_mask) >= config.min_baseline_points
            [row.DisplayVerticalDb, row.DisplayHorizontalDb] = shape_metrics( ...
                map(freq_mask,display_mask),times(display_mask),freqs(freq_mask),config);
            visible = map(:,display_mask);
            finite = isfinite(visible);
            row.DisplayClippedFraction = sum(abs(visible(finite)) >= config.display_limit_db) / ...
                max(1,sum(finite(:))); % 記録のみ。判定には一切使用しない
        end
        values = [row.BaselineVerticalDb row.BaselineHorizontalDb row.BaselineRoughnessDbPer10Ms];
        review_limits = [config.vertical_review_db config.horizontal_review_db config.roughness_review_db_per_10ms];
        severe_limits = [config.vertical_severe_db config.horizontal_severe_db config.roughness_severe_db_per_10ms];
        labels = ["BaselineVertical","BaselineHorizontal","BaselineRoughness"];
        row.ReviewMetricCount = sum(values >= review_limits);
        row.SevereMetricCount = sum(values >= severe_limits);
        if any(~isfinite(values))
            row = raise_candidate(row,"Review");
            row = add_reason(row,"ShapeMetricUnavailable");
        end
        for k = find(values >= review_limits)
            row = add_reason(row,labels(k) + "_Provisional");
        end
        if row.SevereMetricCount >= config.severe_metrics_for_exclude
            row = raise_candidate(row,"ExcludeCandidate");
        elseif row.ReviewMetricCount > 0
            row = raise_candidate(row,"Review");
        end
        if all(isfinite(map(:))) && max(map(:)) - min(map(:)) < 1e-10
            row = raise_candidate(row,"Review");
            row = add_reason(row,"FlatERSPMap_ConfirmEEG");
        end
    catch err
        row = raise_candidate(row,"ExcludeCandidate");
        row = add_reason(row,"ReadOrSchemaError: " + string(err.message));
    end
end

function value = scalar_text(value)
    value = string(value);
    if ~isscalar(value) || ismissing(value)
        value = "";
    end
end

function data = load_available(path,names)
    info = whos('-file',char(path));
    names = intersect(names,{info.name},'stable');
    if isempty(names)
        data = struct();
    else
        data = load(char(path),names{:});
    end
end

function condition = normalize_condition(text_value)
    text_value = scalar_text(text_value);
    if contains(text_value,"NoGreen",'IgnoreCase',true) || contains(text_value,"緑なし")
        condition = "NoGreen";
        return;
    end
    token = regexp(char(text_value), '(?:^|_)(0|80|160)\s*(?:[Oo]|Hz)(?:_|$)', 'tokens','once','ignorecase');
    if isempty(token)
        condition = "Unknown";
    else
        condition = string(token{1}) + " Hz";
    end
end

function [map,times,freqs] = checked_map(data)
    if ~all(isfield(data,{'ersp_mean','times','freqs'}))
        error('erspQC:MissingVariables','ersp_mean/times/freqsがありません。');
    end
    map = data.ersp_mean; times = data.times; freqs = data.freqs;
    if ~isnumeric(map) || ~ismatrix(map) || ~isreal(map) || ...
            ~isnumeric(times) || ~isnumeric(freqs) || ~isvector(times) || ~isvector(freqs) || ...
            ~isreal(times) || ~isreal(freqs)
        error('erspQC:InvalidMap','ERSPは実数行列、軸は実数ベクトルである必要があります。');
    end
    map = double(map); times = double(times(:)'); freqs = double(freqs(:));
    if ~isequal(size(map),[numel(freqs),numel(times)]) || ...
            numel(times) < 3 || numel(freqs) < 3 || ...
            any(~isfinite(times)) || any(~isfinite(freqs)) || ...
            any(diff(times) <= 0) || any(diff(freqs) <= 0) || any(freqs <= 0)
        error('erspQC:InvalidAxes','ERSP寸法と昇順の時間・正の周波数軸が一致しません。');
    end
end

function [vertical,horizontal,roughness] = shape_metrics(map,times,freqs,config)
% 縦筋: 中心列と左右10-30ms中央値の平均との差を、周波数方向で中央値。
%       その時間方向最大値（dB）。左右の窓が揃う場所だけを測る。
% 横筋: 中心行と上下3-10Hz中央値の平均との差を、時間方向で中央値。
%       その周波数方向最大値（dB）。両窓が揃う場所だけを測る。
% 粗さ: 隣接列の|dB差|/時間差を全点で中央値にし10ms当たりに換算。
%       平均ERSPの不規則性であって、EEGのノイズ量の直接推定ではない。
    temporal = nan(1,numel(times));
    for t = 1:numel(times)
        left = times >= times(t)-config.temporal_flank_ms(2) & ...
            times <= times(t)-config.temporal_flank_ms(1);
        right = times >= times(t)+config.temporal_flank_ms(1) & ...
            times <= times(t)+config.temporal_flank_ms(2);
        if any(left) && any(right)
            ref = (median(map(:,left),2,'omitnan') + median(map(:,right),2,'omitnan'))/2;
            temporal(t) = median(abs(map(:,t)-ref),'omitnan');
        end
    end
    spectral = nan(numel(freqs),1);
    for f = 1:numel(freqs)
        lower = freqs >= freqs(f)-config.frequency_flank_hz(2) & ...
            freqs <= freqs(f)-config.frequency_flank_hz(1);
        upper = freqs >= freqs(f)+config.frequency_flank_hz(1) & ...
            freqs <= freqs(f)+config.frequency_flank_hz(2);
        if any(lower) && any(upper)
            ref = (median(map(lower,:),1,'omitnan') + median(map(upper,:),1,'omitnan'))/2;
            spectral(f) = median(abs(map(f,:)-ref),'omitnan');
        end
    end
    vertical = finite_max(temporal);
    horizontal = finite_max(spectral);
    changes = abs(diff(map,1,2))./diff(times)*10;
    roughness = median(changes(:),'omitnan');
end

function value = finite_max(values)
    values = values(isfinite(values));
    if isempty(values), value = NaN; else, value = max(values); end
end

function row = add_reason(row,reason)
    if strlength(row.Reasons) == 0
        row.Reasons = reason;
    else
        row.Reasons = row.Reasons + "; " + reason;
    end
end

function row = set_candidate(row,code)
    row.CandidateCode = code;
    row.CandidateLabel = candidate_label(code);
end

function row = raise_candidate(row,code)
    if candidate_rank(code) > candidate_rank(row.CandidateCode)
        row = set_candidate(row,code);
    end
end

function rank = candidate_rank(code)
    rank = find(["AcceptCandidate","Review","ExcludeCandidate"] == code,1);
    if isempty(rank), rank = 2; end
end

function label = candidate_label(code)
    labels = ["採用候補","要確認","除外候補"];
    label = labels(candidate_rank(code));
end

function cells = flag_duplicates_and_channel_outliers(cells,config)
    keys = cells.ParticipantID + "|" + cells.Condition + "|" + cells.Channel;
    [~,~,groups] = unique(keys);
    for g = unique(groups)'
        indices = find(groups == g);
        if numel(indices) > 1
            for i = indices'
                cells = flag_table_row(cells,i,"ExcludeCandidate","DuplicateCell_DoNotMerge");
            end
        end
    end
    pairs = unique(cells.ParticipantID + "|" + cells.Condition,'stable');
    for p = 1:numel(pairs)
        indices = find((cells.ParticipantID + "|" + cells.Condition) == pairs(p));
        values = cells.BaselineRoughnessDbPer10Ms(indices);
        valid = isfinite(values) & values > 0 & cells.NonfiniteFraction(indices) == 0;
        if sum(valid) < config.min_channels_for_outlier || ...
                numel(unique(cells.Channel(indices(valid)))) < config.min_channels_for_outlier
            continue;
        end
        center = median(values(valid));
        log_values = log10(values(valid));
        log_center = median(log_values);
        scale = max(1.4826*median(abs(log_values-log_center)),config.channel_log_mad_floor);
        for j = find(valid)'
            i = indices(j);
            cells.ChannelRoughnessRatio(i) = values(j)/center;
            cells.ChannelRobustZ(i) = (log10(values(j))-log_center)/scale;
            if cells.ChannelRoughnessRatio(i) >= config.channel_outlier_ratio && ...
                    cells.ChannelRobustZ(i) >= config.channel_outlier_z
                cells = flag_table_row(cells,i,"Review","ChannelRoughnessOutlier_Provisional");
            end
        end
    end
end

function cells = flag_table_row(cells,i,code,reason)
    if candidate_rank(code) > candidate_rank(cells.CandidateCode(i))
        cells.CandidateCode(i) = code;
        cells.CandidateLabel(i) = candidate_label(code);
    end
    if strlength(cells.Reasons(i)) == 0
        cells.Reasons(i) = reason;
    else
        cells.Reasons(i) = cells.Reasons(i) + "; " + reason;
    end
end

function summary = summarize_subjects(subjects,cells,config)
    template = struct('ParticipantID',"",'ExpectedCells',0,'ObservedUniqueCells',0, ...
        'MissingCells',0,'MissingCellList',"",'AcceptCells',0,'ReviewCells',0, ...
        'ExcludeCandidateCells',0,'ExcludeCandidateFraction',0, ...
        'CandidateCode',"Review",'CandidateLabel',"要確認",'Reasons',"", ...
        'FinalDecision',"",'Reviewer',"",'FinalReason',"",'ReviewBoard',"",'RenderError',"");
    rows = repmat(template,height(subjects),1);
    expected = strings(0,1);
    for c = config.conditions
        expected = [expected; c + "|" + config.channels(:)]; %#ok<AGROW>
    end
    for s = 1:height(subjects)
        r = template; r.ParticipantID = subjects.ParticipantID(s);
        selected = cells(cells.ParticipantID == r.ParticipantID,:);
        observed = unique(selected.Condition + "|" + selected.Channel);
        missing = setdiff(expected,observed,'stable');
        r.ExpectedCells = numel(expected);
        r.ObservedUniqueCells = sum(ismember(expected,observed));
        r.MissingCells = numel(missing); r.MissingCellList = strjoin(missing,'; ');
        r.AcceptCells = sum(selected.CandidateCode == "AcceptCandidate");
        r.ReviewCells = sum(selected.CandidateCode == "Review");
        r.ExcludeCandidateCells = sum(selected.CandidateCode == "ExcludeCandidate");
        % 重複セルが候補割合を水増ししないようunique keyで数える。
        bad_keys = unique(selected.Condition(selected.CandidateCode == "ExcludeCandidate") + ...
            "|" + selected.Channel(selected.CandidateCode == "ExcludeCandidate"));
        r.ExcludeCandidateFraction = sum(ismember(expected,bad_keys))/r.ExpectedCells;
        r.CandidateCode = "AcceptCandidate";
        if r.ExcludeCandidateFraction >= config.subject_exclude_cell_fraction
            r.CandidateCode = "ExcludeCandidate";
            r.Reasons = "ManyExcludeCandidateCells_Provisional";
        elseif any(selected.CandidateCode ~= "AcceptCandidate") || r.MissingCells > 0
            r.CandidateCode = "Review";
            r.Reasons = "FlaggedOrMissingCells";
        end
        r.CandidateLabel = candidate_label(r.CandidateCode);
        rows(s) = r;
    end
    summary = struct2table(rows);
end

function output = unique_output_dir(parent)
    if ~isfolder(parent)
        [ok,msg] = mkdir(parent);
        if ~ok, error('erspQC:CannotWrite','%s',msg); end
    end
    base = "run_" + string(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    output = string(fullfile(parent,base));
    suffix = 1;
    while isfolder(output) || isfile(output)
        output = string(fullfile(parent,base + "_" + string(suffix)));
        suffix = suffix + 1;
    end
    [ok,msg] = mkdir(output);
    if ~ok, error('erspQC:CannotWrite','%s',msg); end
end

function render_cell(row,path,config)
    fig = figure('Visible','off','Color','w','Position',[50 50 820 570]);
    cleaner = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig,'Position',[0.10 0.19 0.74 0.68]);
    draw_map(ax,row.SourceMAT,config);
    colormap(fig,jet(256));
    bar = colorbar(ax); bar.Label.String = 'ERSP (dB)';
    title(ax,row.ParticipantID + " / " + row.Condition + " / " + row.Channel + ...
        " / " + row.CandidateLabel,'Interpreter','none');
    xlabel(ax,'Time (ms)'); ylabel(ax,'Frequency (Hz)');
    annotation(fig,'textbox',[0.06 0.01 0.90 0.12], ...
        'String',sprintf('Baseline V=%.3f / H=%.3f dB / rough=%.3f dB/10ms / trials=%g\n%s', ...
        row.BaselineVerticalDb,row.BaselineHorizontalDb,row.BaselineRoughnessDbPer10Ms, ...
        row.TrialCount,char(row.Reasons)),'Interpreter','none','EdgeColor','none','FontSize',9);
    exportgraphics(fig,path,'Resolution',100);
end

function draw_map(ax,path,config)
    try
        data = load_available(path,{'ersp_mean','times','freqs'});
        [map,times,freqs] = checked_map(data);
        use_t = times >= config.display_window_ms(1) & times <= config.display_window_ms(2);
        use_f = freqs >= config.display_frequency_hz(1) & freqs <= config.display_frequency_hz(2);
        if sum(use_t) < 2 || sum(use_f) < 2
            error('erspQC:DisplayCoverage','表示範囲のデータが不足しています。');
        end
        times = times(use_t); freqs = freqs(use_f); map = map(use_f,use_t);
        % imagescの等間隔仮定を避け、実際の非等間隔軸をセル辺で表示。
        % flat shadingで補間しない。末尾行/列も落とさない。
        te = cell_edges(times); fe = cell_edges(freqs);
        padded = map([1:end end],[1:end end]);
        surface(ax,te,fe,zeros(numel(fe),numel(te)),padded, ...
            'FaceColor','flat','EdgeColor','none');
        view(ax,2); set(ax,'YDir','normal');
        xlim(ax,config.display_window_ms); ylim(ax,config.display_frequency_hz);
        caxis(ax,[-config.display_limit_db config.display_limit_db]);
        xline(ax,0,'k--'); box(ax,'on');
    catch err
        axis(ax,'off');
        text(ax,0.05,0.5,"表示不能: " + string(err.message), ...
            'Units','normalized','Interpreter','none','FontSize',8);
    end
end

function edges = cell_edges(values)
    values = values(:)';
    mid = (values(1:end-1)+values(2:end))/2;
    edges = [values(1)-(values(2)-values(1))/2, mid, values(end)+(values(end)-values(end-1))/2];
end

function render_subject_board(id,cells,path,config)
    nc = numel(config.conditions); nch = numel(config.channels);
    fig = figure('Visible','off','Color','w','Position',[10 10 max(900,330*nc) max(650,210*nch)]);
    cleaner = onCleanup(@() close(fig)); %#ok<NASGU>
    layout = tiledlayout(fig,nch,nc,'TileSpacing','compact','Padding','compact');
    title(layout,id + " | 全条件・chの参考一覧（判定は候補） | " + ...
        string(sprintf('color: +/-%g dB',config.display_limit_db)), 'Interpreter','none');
    selected = cells(cells.ParticipantID == id,:);
    for ch = 1:nch
        for c = 1:nc
            ax = nexttile(layout);
            rows = find(selected.Condition == config.conditions(c) & selected.Channel == config.channels(ch));
            if numel(rows) == 1
                row = selected(rows,:);
                draw_map(ax,row.SourceMAT,config);
                label = config.conditions(c) + " / " + config.channels(ch) + " / " + row.CandidateLabel;
            else
                axis(ax,'off');
                label = config.conditions(c) + " / " + config.channels(ch);
                text(ax,0.1,0.5,sprintf('Missing/duplicate: %d files',numel(rows)), ...
                    'Units','normalized','Interpreter','none');
            end
            caxis(ax,[-config.display_limit_db config.display_limit_db]);
            title(ax,label,'Interpreter','none','FontSize',9);
            if ch == nch, xlabel(ax,'Time (ms)'); end
            if c == 1, ylabel(ax,'Frequency (Hz)'); end
        end
    end
    colormap(fig,jet(256));
    bar = colorbar(ax); bar.Layout.Tile = 'east'; bar.Label.String = 'ERSP (dB)';
    exportgraphics(fig,path,'Resolution',100);
end

function write_html(path,cells,subjects,config,anova_info)
    fid = fopen(path,'w','n','UTF-8');
    if fid < 0, error('erspQC:CannotWrite','HTMLを書き込めません。'); end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'<!doctype html><html lang="ja"><meta charset="utf-8"><title>ERSP QC prototype</title>');
    fprintf(fid,['<style>body{font-family:Arial,"Meiryo",sans-serif;margin:24px;background:#f5f7fa;color:#172033}' ...
        'table{border-collapse:collapse;background:white;font-size:13px}th,td{padding:8px;border:1px solid #ccd4df}' ...
        'a{color:#1461a8}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px}' ...
        '.card{padding:12px;background:white;border:1px solid #cbd5e1;border-radius:8px}' ...
        '.card img{width:100%%;height:auto}.AcceptCandidate{border-top:6px solid #4b9b72}' ...
        '.Review{border-top:6px solid #d69c22}' ...
        '.ExcludeCandidate{border-top:6px solid #b95a5a}code{overflow-wrap:anywhere}' ...
        'input,select{padding:8px;margin:4px}summary{cursor:pointer}h2{margin-top:30px}</style>']);
    fprintf(fid,'<h1>ERSP品質スクリーニング・試作版</h1>');
    fprintf(fid,['<p><strong>候補のみ。自動除外・ICA再計算・既存ANOVA変更はしていません。</strong><br>' ...
        '採用候補＝本指標で警告なし。記録品質の保証ではありません。' ...
        '形状の閾値は未検証の仮値で、本来の反応を誤検出する可能性があります。' ...
        'イベント後の反応の大きさ、赤い面積、条件差は判定に使っていません。</p>']);
    fprintf(fid,'<p>画像は共通 ±%g dB、時間 %g〜%g ms、周波数 %g〜%g Hz。表示補間なし。</p>', ...
        config.display_limit_db,config.display_window_ms,config.display_frequency_hz);
    fprintf(fid,'<p><a href="QC_cells.csv">セル別指標CSV</a> / <a href="QC_participants.csv">被験者別候補CSV</a></p>');
    if anova_info.Status=="Complete"
        fprintf(fid,'<p>新規07 ANOVA表（候補判定では自動除外せず、完全データを集計）：');
        for k=1:numel(anova_info.Files)
            [~,name,ext]=fileparts(anova_info.Files(k));
            fprintf(fid,' <a href="%s">%s</a>',html_escape(file_url(anova_info.Files(k))),html_escape(name+ext));
        end
        fprintf(fid,'</p>');
    elseif anova_info.Status=="Failed"
        fprintf(fid,'<p>07 ANOVA表の作成に失敗：%s</p>',html_escape(anova_info.Error));
    end
    fprintf(fid,'<h2>被験者一覧</h2><table><tr><th>ID</th><th>候補</th><th>採用候補セル</th><th>要確認</th><th>除外候補</th><th>欠測</th><th>全条件・ch一覧</th></tr>');
    for s = 1:height(subjects)
        r = subjects(s,:);
        board = "—";
        if strlength(r.ReviewBoard) > 0
            board = "<a href='" + html_escape(relative_url(r.ReviewBoard)) + "'>一覧画像</a>";
        end
        fprintf(fid,'<tr><td>%s</td><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td title="%s">%d</td><td>%s</td></tr>', ...
            html_escape(r.ParticipantID),html_escape(r.CandidateLabel),r.AcceptCells,r.ReviewCells, ...
            r.ExcludeCandidateCells,html_escape(r.MissingCellList),r.MissingCells,board);
    end
    fprintf(fid,'</table><h2>全条件・チャンネルの画像</h2>');
    fprintf(fid,'<p>「全て」は採用候補も含む全セル。「両方」は要確認＋除外候補です。</p>');
    fprintf(fid,['<input id="q" placeholder="被験者・条件・chで検索" oninput="filterCards()">' ...
        '<select id="status" onchange="filterCards()"><option value="All">全て</option>' ...
        '<option value="Flagged">両方</option>' ...
        '<option value="Review">要確認</option><option value="ExcludeCandidate">除外候補</option></select>' ...
        '<span id="count"></span><div class="cards">']);
    card_indices = (1:height(cells))';
    [~,order] = sort(arrayfun(@(i) -candidate_rank(cells.CandidateCode(i)),card_indices));
    card_indices = card_indices(order);
    for i = card_indices'
        r = cells(i,:);
        heading = r.ParticipantID + " / " + r.Condition + " / " + r.Channel;
        fprintf(fid,'<article class="card %s" data-status="%s" data-search="%s"><h3>%s</h3><p>%s — %s</p>', ...
            r.CandidateCode,r.CandidateCode,html_escape(lower(heading)),html_escape(heading), ...
            html_escape(r.CandidateLabel),html_escape(r.FileID));
        if strlength(r.ReviewImage) > 0
            fprintf(fid,'<a href="%s"><img loading="lazy" src="%s" alt="%s"></a>', ...
                html_escape(r.ReviewImage),html_escape(r.ReviewImage),html_escape(heading));
        else
            fprintf(fid,'<p>画像なし: %s</p>',html_escape(r.RenderError));
        end
        fprintf(fid,'<p>試行数 %g / baseline縦筋 %.3f dB / 横筋 %.3f dB / 粗さ %.3f dB/10ms</p><p><code>%s</code></p>', ...
            r.TrialCount,r.BaselineVerticalDb,r.BaselineHorizontalDb,r.BaselineRoughnessDbPer10Ms, ...
            html_escape(r.Reasons));
        fprintf(fid,'<details><summary>元ファイル（表示色範囲は元の設定）</summary><p><code>%s</code></p>',html_escape(r.SourceMAT));
        if isfile(r.SourcePNG)
            source_url = file_url(r.SourcePNG);
            fprintf(fid,'<a href="%s">元ERSP画像</a>',html_escape(source_url));
        end
        fprintf(fid,'</details></article>');
    end
    fprintf(fid,['</div><script>function filterCards(){let q=document.getElementById("q").value.toLowerCase(),' ...
        's=document.getElementById("status").value,n=0;document.querySelectorAll(".card").forEach(c=>{' ...
        'let statusMatches=s==="All"||(s==="Flagged"&&(c.dataset.status==="Review"||' ...
        'c.dataset.status==="ExcludeCandidate"))||c.dataset.status===s;' ...
        'let show=c.dataset.search.includes(q)&&statusMatches;c.hidden=!show;if(show)n++;});' ...
        'document.getElementById("count").textContent=n+" 件";}filterCards();</script>']);
    fprintf(fid,'<h2>指標・閾値の定義</h2><pre>%s</pre>',html_escape(string(evalc('disp(config)'))));
    fprintf(fid,['<p>baseline内の周辺列/行との差を使用。非有限値、少数試行、欠測・重複を別に記録。' ...
        '群間のERSP外れ値や統計の有意性で判定しない。' ...
        'このレポートではPSD・電極抵抗・個々の試行は評価していません。</p></html>']);
end

function value = html_escape(value)
    value = string(value);
    value = replace(value,'&','&amp;'); value = replace(value,'<','&lt;');
    value = replace(value,'>','&gt;'); value = replace(value,'"','&quot;');
    value = replace(value,"'",'&#39;');
end

function url = file_url(path)
    url = string(char(java.io.File(char(path)).toURI().toASCIIString()));
end

function url = relative_url(path)
% 被験者名の空白、#、%、&、日本語もブラウザで正しいファイルを指すようにする。
    parts = split(replace(string(path), "\", "/"), "/");
    for k = 1:numel(parts)
        parts(k) = replace(string(java.net.URLEncoder.encode(char(parts(k)), 'UTF-8')), '+', '%20');
    end
    url = strjoin(parts, '/');
end

function cells = empty_threeway_cells()
    cells = table('Size', [0 9], ...
        'VariableTypes', {'string', 'string', 'string', 'string', ...
            'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'ParticipantID', 'Condition', 'Band', 'Channel', ...
            'ERSPdB', 'TimePointCount', 'MinimumTrialCount', ...
            'WindowStartMs', 'WindowEndMs'});
end

function config = threeway_factor_config(config)
    config.anova_conditions=config.conditions;
    config.band_names=config.anova_band_names;
    config.analysis_channels=config.channels;
    [~,index]=ismember(config.conditions,["0 Hz","80 Hz","160 Hz","NoGreen"]);
    tags=["C0Hz","C80Hz","C160Hz","NoGreen"];
    config.anova_condition_tags=tags(index);
end

function info = export_band_anova(subjects,quality,output_dir,config)
% 06 MAT優先。06がない条件だけ04の保存済みERSPから同じ帯域平均を復元する。
% ICA/newtimef/dB変換を再実行せず、不明な有効試行数は推測しない。
    config=threeway_factor_config(config);
    cells=empty_threeway_cells();
    sources=table('Size',[0 8], ...
        'VariableTypes',{'string','string','string','string','string','double','double','string'}, ...
        'VariableNames',{'ParticipantID','Condition','Channel','SourcePath','Status', ...
            'ICASeed','BrainSuggestionThreshold','ICASelectionMode'});
    for s=1:height(subjects)
        for ch=config.channels
            folder=fullfile(subjects.SubjectFolder(s),config.band_folder_name,ch);
            mat_files=dir(fullfile(folder,'*_共通ICA_帯域ERSP.mat'));
            csv_files=dir(fullfile(folder,'*_共通ICA_帯域ERSP.csv'));
            paths=string(fullfile({mat_files.folder},{mat_files.name}));
            % 同じstemのMAT/CSVを二重計上しない。別stemの同一セルは重複として残す。
            for f=1:numel(csv_files)
                path=string(fullfile(csv_files(f).folder,csv_files(f).name));
                if ~ismember(replace(path,'.csv','.mat'),paths), paths(end+1)=path; end %#ok<AGROW>
            end
            paths=sort(paths);
            for path=paths
                [rows,meta]=read_band_cells(path,subjects.ParticipantID(s),ch,config);
                if isempty(meta), continue; end
                if ~isempty(rows) && ~ismember(rows.Condition(1),config.conditions), continue; end
                cells=[cells;rows]; %#ok<AGROW>
                sources=[sources;meta]; %#ok<AGROW>
            end
            % 06がない条件だけ04を使用。読取不能な06を黙って04へ置き換えない。
            primary_folder=fullfile(subjects.SubjectFolder(s),config.ersp_folder_name,ch);
            primary_files=dir(fullfile(primary_folder,'*_共通ICA_newtimef_ERSP.mat'));
            band_sources=sources(sources.ParticipantID==subjects.ParticipantID(s) & sources.Channel==ch,:);
            for f=1:numel(primary_files)
                path=string(fullfile(primary_files(f).folder,primary_files(f).name));
                [rows,meta]=read_ersp_band_cells(path,subjects.ParticipantID(s),ch,config);
                if isempty(meta), continue; end
                already=band_sources.Condition==meta.Condition;
                if any(already), continue; end
                cells=[cells;rows]; sources=[sources;meta]; %#ok<AGROW>
            end
        end
    end
    if isempty(sources)
        error('erspQC:NoBandFiles','選択した版の04/06解析ファイルが0件です。入力範囲を確認してください。');
    end
    info=save_anova_tables(cells,subjects.ParticipantID,output_dir,config, ...
        "ExistingBandOutputs_04FallbackWhen06Missing",sources,quality);
end

function [rows,meta] = read_band_cells(path,id,channel,config)
    rows=empty_threeway_cells();
    condition="Unknown";
    seed=NaN; brain=NaN; selection="Unknown";
    try
        [~,base,ext]=fileparts(path);
        condition=normalize_condition(string(base));
        if strcmpi(ext,'.mat')
            data=load_available(path,{'band_mean','band_n','times','condition_label','analysis_config'});
            if isfield(data,'condition_label'), condition=normalize_condition(scalar_text(data.condition_label)); end
            if condition=="Unknown", condition=normalize_condition(string(base)); end
            if ~ismember(condition,config.conditions)
                if condition=="Unknown", error('erspQC:UnknownBandCondition','06の条件を判定できません。'); end
                meta=table(); return;
            end
            if ~all(isfield(data,{'band_mean','band_n','times','analysis_config'})) || ...
                    ~isfield(data.analysis_config,'band_names')
                error('erspQC:BandSchema','06 MATの帯域・時間・N・帯域名が不足しています。');
            end
            names=string(data.analysis_config.band_names(:));
            times=double(data.times(:)'); values=double(data.band_mean); counts=double(data.band_n);
            if ~isreal(values) || ~isreal(counts) || ...
                    ~isequal(size(values),[numel(names),numel(times)]) || ~isequal(size(counts),size(values)) || ...
                    numel(unique(names))~=numel(names)
                error('erspQC:BandShape','06 MATの帯域行・時間列・Nのサイズが不正です。');
            end
            cfg=data.analysis_config;
            if isfield(cfg,'ica_seed') && isnumeric(cfg.ica_seed) && isscalar(cfg.ica_seed), seed=double(cfg.ica_seed); end
            if isfield(cfg,'brain_keep_threshold') && isnumeric(cfg.brain_keep_threshold) && isscalar(cfg.brain_keep_threshold)
                brain=double(cfg.brain_keep_threshold);
            end
            if isfield(cfg,'ica_selection_mode'), selection=scalar_text(cfg.ica_selection_mode); end
            for band=config.band_names
                row=blank_band_cell(id,condition,band,channel,config);
                index=find(names==band);
                if numel(index)==1
                    row=summarize_band_vector(row,times,values(index,:),counts(index,:),config);
                end
                rows=[rows;row]; %#ok<AGROW>
            end
            status="BandMAT";
        else
            source=readtable(path,'TextType','string','Encoding','UTF-8');
            required={'Condition','Channel','Band','TimeMs','MeanERSPdB','TrialCount'};
            if ~all(ismember(required,source.Properties.VariableNames))
                error('erspQC:BandCSVSchema','06 CSVの必須列が不足しています。');
            end
            labels=unique(string(source.Condition));
            labels=arrayfun(@normalize_condition,labels);
            if numel(labels)~=1, error('erspQC:MixedBandConditions','06 CSVに複数条件があります。'); end
            condition=labels(1);
            if condition=="Unknown", condition=normalize_condition(string(base)); end
            if ~ismember(condition,config.conditions)
                if condition=="Unknown", error('erspQC:UnknownBandCondition','06 CSVの条件を判定できません。'); end
                meta=table(); return;
            end
            if any(string(source.Channel)~=channel), error('erspQC:BandCSVChannel','06 CSVのchとフォルダ名が一致しません。'); end
            for band=config.band_names
                selected=source(string(source.Band)==band,:);
                [~,order]=sort(selected.TimeMs); selected=selected(order,:);
                row=blank_band_cell(id,condition,band,channel,config);
                if ~isempty(selected)
                    row=summarize_band_vector(row,selected.TimeMs',selected.MeanERSPdB',selected.TrialCount',config);
                end
                rows=[rows;row]; %#ok<AGROW>
            end
            status="BandCSV_ICASettingsUnknown";
        end
    catch err
        % 不正・読取不能は採用せず、欠測/invalidとして理由と入力参照を記録する。
        rows=empty_threeway_cells();
        if ismember(condition,config.conditions)
            for band=config.band_names, rows=[rows;blank_band_cell(id,condition,band,channel,config)]; end %#ok<AGROW>
        end
        status="ReadOrSchemaError: "+string(err.message);
    end
    meta=table(string(id),condition,string(channel),string(path),status,seed,brain,selection, ...
        'VariableNames',{'ParticipantID','Condition','Channel','SourcePath','Status', ...
            'ICASeed','BrainSuggestionThreshold','ICASelectionMode'});
end

function [rows,meta] = read_ersp_band_cells(path,id,channel,config)
% align_band_means_with_newtimefと同じ周波数方向dB平均。
% 04のNだけでは試行の周波数間の欠測対応を一般には復元できない。
% 指定窓の全周波数でN==trial_countのときに限り、band_nが厳密に求まる。
    rows=empty_threeway_cells(); condition="Unknown";
    seed=NaN; brain=NaN; selection="Unknown";
    try
        [~,base]=fileparts(path); condition=normalize_condition(string(base));
        data=load_available(path,{'ersp_mean','times','freqs','ersp_n','trial_count','condition_label','analysis_config'});
        if isfield(data,'condition_label')
            label=normalize_condition(scalar_text(data.condition_label));
            if label~="Unknown", condition=label; end
        end
        if ~ismember(condition,config.conditions)
            if condition=="Unknown", error('erspQC:UnknownERSPCondition','04の条件が不明です。'); end
            meta=table(); return;
        end
        [map,times,freqs]=checked_map(data);
        if ~all(isfield(data,{'ersp_n','trial_count','analysis_config'})) || ...
                ~all(isfield(data.analysis_config,{'band_names','band_ranges'}))
            error('erspQC:ERSPBandSchema','04から帯域平均・試行数を復元する情報が不足しています。');
        end
        cfg=data.analysis_config; names=string(cfg.band_names(:)); ranges=double(cfg.band_ranges);
        if size(ranges,1)~=numel(names) || size(ranges,2)~=2 || numel(unique(names))~=numel(names) || ...
                any(~isfinite(ranges(:))) || any(ranges(:,2)<ranges(:,1))
            error('erspQC:ERSPBandRanges','04に保存された帯域定義が不正です。');
        end
        n=data.ersp_n; trials=data.trial_count;
        if ~isnumeric(n) || ~isreal(n) || ~isequal(size(n),size(map)) || ...
                ~isnumeric(trials) || ~isscalar(trials) || ~isfinite(trials) || trials<1 || trials~=fix(trials)
            error('erspQC:ERSPBandCounts','04に保存された試行数が不正です。');
        end
        if isfield(cfg,'ica_seed') && isnumeric(cfg.ica_seed) && isscalar(cfg.ica_seed), seed=double(cfg.ica_seed); end
        if isfield(cfg,'brain_keep_threshold') && isnumeric(cfg.brain_keep_threshold) && isscalar(cfg.brain_keep_threshold)
            brain=double(cfg.brain_keep_threshold);
        end
        if isfield(cfg,'ica_selection_mode'), selection=scalar_text(cfg.ica_selection_mode); end
        window=times>=config.anova_time_window_ms(1) & times<=config.anova_time_window_ms(2);
        for band=config.band_names
            row=blank_band_cell(id,condition,band,channel,config);
            index=find(names==band);
            if numel(index)==1
                mask=freqs>=ranges(index,1) & freqs<=ranges(index,2);
                if any(mask)
                    if any(n(mask,window)~=trials,'all')
                        error('erspQC:AmbiguousBandTrialCount','04だけでは欠測を含む帯域の有効試行数を確定できません。06が必要です。');
                    end
                    means=mean(map(mask,:),1,'omitnan');
                    row=summarize_band_vector(row,times,means,repmat(double(trials),1,numel(times)),config);
                end
            end
            rows=[rows;row]; %#ok<AGROW>
        end
        status="ERSPMAT_DerivedBandMeans_FullTrialCounts";
    catch err
        rows=empty_threeway_cells();
        if ismember(condition,config.conditions)
            for band=config.band_names, rows=[rows;blank_band_cell(id,condition,band,channel,config)]; end %#ok<AGROW>
        end
        status="ReadOrSchemaError: "+string(err.message);
    end
    meta=table(string(id),condition,string(channel),string(path),status,seed,brain,selection, ...
        'VariableNames',{'ParticipantID','Condition','Channel','SourcePath','Status', ...
            'ICASeed','BrainSuggestionThreshold','ICASelectionMode'});
end

function row = blank_band_cell(id,condition,band,channel,config)
    template=empty_threeway_cells();
    row=table(string(id),string(condition),string(band),string(channel),NaN,0,NaN, ...
        config.anova_time_window_ms(1),config.anova_time_window_ms(2), ...
        'VariableNames',template.Properties.VariableNames);
end

function row = summarize_band_vector(row,times,values,counts,config)
    if ~isnumeric(times) || ~isreal(times) || numel(times)<2 || ...
            any(~isfinite(times)) || any(diff(times)<=0) || ...
            times(1)>config.anova_time_window_ms(1) || times(end)<config.anova_time_window_ms(2)
        return;
    end
    selected=times>=config.anova_time_window_ms(1) & times<=config.anova_time_window_ms(2);
    row.TimePointCount=sum(selected);
    if numel(values)~=numel(times) || numel(counts)~=numel(times), return; end
    values=values(selected); counts=counts(selected);
    if ~isempty(values) && all(isfinite(values)), row.ERSPdB=mean(values); end
    if ~isempty(counts) && all(isfinite(counts)) && all(counts>=1) && all(counts==fix(counts))
        row.MinimumTrialCount=min(counts);
    end
end

function info = save_anova_tables(cells,ids,output_dir,config,origin,sources,quality)
    [wide,design,inclusion,complete_long]=build_threeway_tables(cells,ids,config);
    inclusion.DataOrigin=repmat(string(origin),height(inclusion),1);
    inclusion.ICASeed=nan(height(inclusion),1);
    inclusion.BrainKeepThreshold=nan(height(inclusion),1);
    inclusion.BrainSuggestionThreshold=nan(height(inclusion),1);
    inclusion.ICASelectionMode=repmat("Unknown",height(inclusion),1);
    inclusion.ERSPQualityCandidate=repmat("NotEvaluated",height(inclusion),1);
    for s=1:height(inclusion)
        id=inclusion.ParticipantID(s);
        if ~isempty(sources)
            selected=sources(sources.ParticipantID==id,:);
            if ~isempty(selected)
                if all(isfinite(selected.ICASeed)) && numel(unique(selected.ICASeed))==1
                    inclusion.ICASeed(s)=selected.ICASeed(1);
                end
                if all(isfinite(selected.BrainSuggestionThreshold)) && numel(unique(selected.BrainSuggestionThreshold))==1
                    inclusion.BrainSuggestionThreshold(s)=selected.BrainSuggestionThreshold(1);
                end
                if numel(unique(selected.ICASelectionMode))==1
                    inclusion.ICASelectionMode(s)=selected.ICASelectionMode(1);
                    if selected.ICASelectionMode(1)=="brain"
                        inclusion.BrainKeepThreshold(s)=inclusion.BrainSuggestionThreshold(s);
                    end
                end
            end
        end
        if ~isempty(quality)
            index=find(quality.ParticipantID==id);
            if numel(index)==1, inclusion.ERSPQualityCandidate(s)=quality.CandidateCode(index); end
        end
    end
    files=string(fullfile(output_dir,["ANOVA_all_cells_QC.csv";"ANOVA_threeway_wide.csv"; ...
        "ANOVA_threeway_factor_design.csv";"ANOVA_threeway_inclusion_QC.csv";"ANOVA_threeway_long_complete.csv"]));
    outputs={cells,wide,design,inclusion,complete_long};
    for k=1:numel(files), writetable(outputs{k},files(k),'Encoding','UTF-8'); end
    info=struct('Status',"Complete",'OutputDir',string(output_dir),'Files',files, ...
        'IncludedParticipants',height(wide),'MeasurementColumns',height(design), ...
        'SourceFiles',sources,'TimeWindowMs',config.anova_time_window_ms, ...
        'QualityFiltering',"None: only completeness/finite/count/window requirements");
    fprintf('07 三要因wide: %d名 / %d列（%d条件×%d帯域×%dch）\n%s\n', ...
        height(wide),height(design),numel(config.conditions),numel(config.band_names),numel(config.channels),output_dir);
    if height(wide)<2
        warning('erspQC:TooFewANOVAParticipants','完全な被験者が2名未満です。群の分散分析はできません。');
    end
end

function [output_dir,info] = export_existing_anova(input_root,output_root,config)
    path=fullfile(input_root,'ANOVA_all_cells_QC.csv');
    if ~isfile(path), path=fullfile(input_root,config.anova_folder_name,'ANOVA_all_cells_QC.csv'); end
    if ~isfile(path), error('erspQC:MissingANOVACells','ANOVA_all_cells_QC.csvのあるフォルダを指定してください。'); end
    source=readtable(path,'TextType','string','Encoding','UTF-8');
    template=empty_threeway_cells(); required=template.Properties.VariableNames;
    if ~all(ismember(required,source.Properties.VariableNames)), error('erspQC:ANOVACellSchema','既存cells表の必須列が不足しています。'); end
    cells=source(:,required);
    cells.Condition=arrayfun(@normalize_condition,cells.Condition);
    cells=cells(ismember(cells.Condition,config.conditions) & ismember(cells.Band,config.anova_band_names) & ...
        ismember(cells.Channel,config.channels),:);
    if isempty(cells), error('erspQC:EmptyANOVACells','指定の因子に対応する既存cells表が空です。'); end
    if strlength(config.anova_output_root)>0, parent=config.anova_output_root;
    elseif strlength(output_root)>0, parent=fullfile(output_root,config.anova_folder_name);
    else, parent=string(fileparts(path)); end
    output_dir=unique_output_dir(parent);
    config=threeway_factor_config(config);
    info=save_anova_tables(cells,unique(cells.ParticipantID,'stable'),output_dir,config, ...
        "ExistingANOVACSV_UnverifiedICASettings",table(),table());
end

function [wide, design, inclusion, complete_long] = ...
        build_threeway_tables(cells, participant_ids, config)
    conditions = config.anova_conditions(:)';
    bands = config.band_names(:)';
    channels = config.analysis_channels(:)';
    tags = config.anova_condition_tags(:)';
    if numel(tags) ~= numel(conditions) || ...
            numel(unique(conditions)) ~= numel(conditions) || ...
            numel(unique(bands)) ~= numel(bands) || ...
            numel(unique(lower(channels))) ~= numel(channels)
        error('erspQC:InvalidANOVAFactors', ...
            '三要因の水準または条件タグが不正・重複しています。');
    end
    validateattributes(config.anova_time_window_ms, {'numeric'}, ...
        {'real', 'finite', 'numel', 2, 'increasing'});
    if any(~ismember(cells.Condition, conditions))
        error('erspQC:UnexpectedANOVACondition', ...
            '三要因表に未設定の条件があります: %s', ...
            strjoin(unique(cells.Condition(~ismember(cells.Condition, ...
                conditions))), ', '));
    end
    participant_ids = string(participant_ids(:));
    if numel(unique(participant_ids)) ~= numel(participant_ids)
        error('erspQC:DuplicateParticipantID', ...
            '三要因用の被験者IDが重複しています。');
    end
    n_factors = numel(conditions) * numel(bands) * numel(channels);
    design = table('Size', [n_factors 8], ...
        'VariableTypes', {'string', 'string', 'string', 'string', ...
            'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Variable', 'Condition', 'Band', 'Channel', ...
            'ConditionOrder', 'BandOrder', 'ChannelOrder', 'MeasurementOrder'});
    position = 0;
    for condition_index = 1:numel(conditions)
        for band_index = 1:numel(bands)
            for channel_index = 1:numel(channels)
                position = position + 1;
                name = tags(condition_index) + "_" + bands(band_index) + ...
                    "_" + channels(channel_index);
                if ~isvarname(char(name)) || strlength(name) > 64
                    error('erspQC:InvalidSPSSVariable', ...
                        'SPSS用の変数名が不正です: %s', name);
                end
                design(position, :) = {name, conditions(condition_index), ...
                    bands(band_index), channels(channel_index), ...
                    condition_index, band_index, channel_index, position};
            end
        end
    end
    if numel(unique(design.Variable)) ~= n_factors
        error('erspQC:DuplicateSPSSVariable', ...
            'SPSS用の測定変数名が重複しています。');
    end
    wide = array2table(zeros(0, n_factors), ...
        'VariableNames', cellstr(design.Variable));
    wide = addvars(wide, strings(0, 1), 'Before', 1, ...
        'NewVariableNames', 'ParticipantID');
    inclusion = table('Size', [0 5], ...
        'VariableTypes', {'string', 'logical', 'double', 'double', 'string'}, ...
        'VariableNames', {'ParticipantID', 'Included', 'ObservedCells', ...
            'ExpectedCells', 'MissingOrInvalidCells'});
    complete_long = cells([], :);
    for participant_index = 1:numel(participant_ids)
        id = participant_ids(participant_index);
        selected = cells(cells.ParticipantID == id & ...
            ismember(cells.Band, bands) & ...
            ismember(cells.Channel, channels), :);
        missing = strings(0, 1);
        values = nan(1, n_factors);
        selected_order = zeros(n_factors, 1);
        expected_time_points = [];
        for factor_index = 1:n_factors
            match = find(selected.Condition == design.Condition(factor_index) & ...
                selected.Band == design.Band(factor_index) & ...
                selected.Channel == design.Channel(factor_index));
            if numel(match) ~= 1
                missing(end + 1, 1) = design.Variable(factor_index) + ...
                    " (count=" + numel(match) + ")"; %#ok<AGROW>
                continue;
            end
            row = selected(match, :);
            if isempty(expected_time_points) && ...
                    isfinite(row.TimePointCount) && row.TimePointCount >= 1
                expected_time_points = row.TimePointCount;
            end
            valid = isfinite(row.ERSPdB) && ...
                isfinite(row.MinimumTrialCount) && row.MinimumTrialCount >= 1 && ...
                isfinite(row.TimePointCount) && row.TimePointCount >= 1 && ...
                row.TimePointCount == fix(row.TimePointCount) && ...
                ~isempty(expected_time_points) && ...
                row.TimePointCount == expected_time_points && ...
                row.WindowStartMs == config.anova_time_window_ms(1) && ...
                row.WindowEndMs == config.anova_time_window_ms(2);
            if ~valid
                missing(end + 1, 1) = design.Variable(factor_index) + ...
                    " (invalid value/count/window)"; %#ok<AGROW>
                continue;
            end
            values(factor_index) = row.ERSPdB;
            selected_order(factor_index) = match;
        end
        included = isempty(missing);
        inclusion(end + 1, :) = {id, included, height(selected), ...
            n_factors, strjoin(missing, '; ')}; %#ok<AGROW>
        if included
            row = array2table(values, ...
                'VariableNames', cellstr(design.Variable));
            row = addvars(row, id, 'Before', 1, ...
                'NewVariableNames', 'ParticipantID');
            wide = [wide; row]; %#ok<AGROW>
            complete_long = [complete_long; selected(selected_order, :)]; %#ok<AGROW>
        end
    end
end
