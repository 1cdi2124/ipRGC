function input_dirs = pae_select_subject_input_dirs( ...
        selected_root, preferred_stage, extension, fallback_stage)
%PAE_SELECT_SUBJECT_INPUT_DIRS Resolve one subject or a direct-child cohort.
% Select an input stage folder, its subject folder, or a folder containing
% subject folders. The returned input folders are sorted by subject name.
% Only direct child subjects are considered, so nested experiments cannot
% be mixed accidentally. A selected stage folder always takes precedence.

    if nargin < 4
        fallback_stage = "";
    end
    selected_root = string(selected_root);
    preferred_stage = string(preferred_stage);
    fallback_stage = string(fallback_stage);
    extension = string(extension);
    if ~isscalar(selected_root) || ~isfolder(selected_root)
        error('pae_select_subject_input_dirs:MissingFolder', ...
            '選択したフォルダが見つかりません: %s', selected_root);
    end
    stage_names = preferred_stage;
    if strlength(fallback_stage) > 0
        stage_names(end + 1) = fallback_stage;
    end

    [~, selected_name] = fileparts(selected_root);
    if ismember(string(selected_name), stage_names)
        require_input_files(selected_root, extension);
        input_dirs = selected_root;
        return;
    end

    stage_dir = find_subject_stage(selected_root, stage_names);
    if strlength(stage_dir) > 0
        require_input_files(stage_dir, extension);
        input_dirs = stage_dir;
        return;
    end

    % Keep the previous behavior for an explicitly selected folder that
    % contains input files directly but does not use the standard name.
    if has_input_files(selected_root, extension)
        input_dirs = selected_root;
        return;
    end

    entries = dir(selected_root);
    entries = entries([entries.isdir]);
    names = string({entries.name});
    names = sort(names(~ismember(names, [".", ".."])));
    input_dirs = strings(0, 1);
    for index = 1:numel(names)
        subject_dir = string(fullfile(selected_root, names(index)));
        stage_dir = find_subject_stage(subject_dir, stage_names);
        if strlength(stage_dir) == 0
            continue;
        end
        require_input_files(stage_dir, extension);
        input_dirs(end + 1, 1) = stage_dir; %#ok<AGROW>
    end
    if isempty(input_dirs)
        error('pae_select_subject_input_dirs:NoSubjects', ...
            '直下に%sを持つ被験者フォルダがありません: %s', ...
            preferred_stage, selected_root);
    end
end

function stage_dir = find_subject_stage(subject_dir, stage_names)
    stage_dir = "";
    for index = 1:numel(stage_names)
        candidate = string(fullfile(subject_dir, stage_names(index)));
        if isfolder(candidate)
            stage_dir = candidate;
            return;
        end
    end
end

function require_input_files(folder, extension)
    if ~has_input_files(folder, extension)
        error('pae_select_subject_input_dirs:NoFiles', ...
            '入力ファイル（*%s）がありません: %s', extension, folder);
    end
end

function found = has_input_files(folder, extension)
    files = dir(fullfile(folder, "*" + extension));
    found = any(~[files.isdir]);
end
