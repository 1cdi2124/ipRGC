function added_paths = setup_analysis_paths(category)
%SETUP_ANALYSIS_PATHS Add the current analysis folder and one optional category.
% Usage: setup_analysis_paths
%        setup_analysis_paths("10_従来方式_CWT")
% Adds paths for this session only; does not run analysis, remove paths, or savepath.
% Backups are deliberately not supported. Do not use genpath on this repository.
    root_dir = fileparts(mfilename('fullpath'));
    allowed = ["10_従来方式_CWT", "20_旧版_共通ICA", ...
        "30_旧ANOVA_二要因", "40_イベント作成", "50_検討用"];
    added_paths = string(root_dir);
    if nargin >= 1 && ~isempty(category)
        category = string(category);
        if ~isscalar(category) || ~ismember(category, allowed)
            error('setup_analysis_paths:InvalidCategory', ...
                'Choose one of: %s', strjoin(allowed, ', '));
        end
        selected_dir = fullfile(root_dir, char(category));
        if ~isfolder(selected_dir)
            error('setup_analysis_paths:MissingFolder', ...
                'Folder does not exist: %s', selected_dir);
        end
        addpath(selected_dir, '-end');
        added_paths(end + 1) = string(selected_dir);
    end
    addpath(root_dir, '-begin');
    fprintf('Analysis path: %s\n', strjoin(added_paths, newline));
end
