function manifest = run_single_ic_ersp_analysis(input_root, analysis_channels, ...
        ica_seed, ersp_max_db, ic_numbers, options)
%RUN_SINGLE_IC_ERSP_ANALYSIS 各ICを1つだけ採用したERSP一覧を出す補助プログラム。
% 通常解析・品質判定・ANOVAは行わない。元SET/04～07/既存QCは変更しない。
% 1 input_root: 03_手動サッケード/被験者/被験者一覧。既定: uigetdir。
% 2 analysis_channels: 出力電極。既定: config.analysis_channels="Oz"。
%   []でOz。複数電極の例: ["O1","Oz","O2"]。ICA学習の8電極は変えない。
% 3 ica_seed: 既定config.ica_seed=3。0～4294967295の整数。
% 4 ersp_max_db: 既定config.ersp_color_limit_db=1.5（画像±1.5dB）。保存値は変更しない。
% 5 ic_numbers: 既定[]=全IC。[2 6]ならIC2単独とIC6単独。2+6同時採用ではない。
%   Brain確率では絞らず、指定した各ICを1つずつ投影する。
% 6 options: default_configの項目を上書きするscalar struct。通常は省略。
%   config.output_parent: 個人選択時は被験者フォルダ/08_ERSP_QCへ自動保存。
%   03_手動サッケードを選んでも、その親の被験者フォルダ内に08を作る。
%   一覧選択時の既定は設定欄の2610/01_CWT/08_ERSP_QC（従来どおり）。
%   options.output_parentを明示した場合は、その指定先を優先する。
%   config.ica_model_file="": 全条件SETを連結して共通ICAを1回学習する。
%   保存済み共通ICA分解MATを指定すれば、その分解のIC番号で出力（1被験者限定）。
%   再学習時、seedが同じでも入力・処理が違えば旧IC番号とは一致しない。
%
% 出力: 08_ERSP_QC/IC_run_日時/images/被験者名_IC01.png, IC02.png, ...
% images: ICごとに全条件×指定電極を表示する個別画像。
% boards/被験者名.png: 出力した全IC×全条件×指定電極をまとめた一覧1枚。
% images/boardsとも、各ICのERSPの左隣に頭皮分布（ICA電極の投影係数）を表示。
% 頭皮分布はIC内の最大絶対値で規格化した±1表示で、ERSPのdBとは別。
% options.ersp_display_max_hz（既定200）で画像の周波数上限だけ変更可能。
% ICの極性・スケールは任意。頭皮分布の赤色はERSP増加を意味しない。
% 座標がない電極のみ標準10-5座標を表示用に補完。解析データは変更しない。
% 失敗/欠測は明示。毎回新規runを作る。既存runの画像は移動・変更しない。
% subjects/被験者名/ にERSPのMAT・使用ICA・設定・処理記録を保持する。
% PSD集計、SET保存、QC分類、HTML、帯域集計、ANOVAは追加しない。
%
% run_single_ic_ersp_analysis
% run_single_ic_ersp_analysis("", "Oz", 3, 1.5, [])
% run_single_ic_ersp_analysis("", ["O1","Oz","O2"], 3, 1.5, 6)
% opt=struct('ica_model_file',"C:\...\ICA_QC\共通ICA_分解と除去判定.mat");
% run_single_ic_ersp_analysis("C:\...\被験者", "Oz", [], 1.5, 6, opt)
%
% CWT_20261006と同じ前処理・boundary保持・共通ICA・newtimefを使用。
% 新規ICA: 0.2～220Hz/平均基準 → 学習1～100Hz/250Hz、上限rank7。
% A(:,k)*W(k,:)*EEGの単独IC投影のみ採用。他IC・残差は足さない。
% Saccade[-4,4]秒、時間baseline[-250,-100]ms、newtimef0.5～200Hz。
% 平均ERSPはnewtimefの返す値を使用。線形補間/独自dB平均はしない。
% 単独ICでは各電極が同一波形の定数倍なので相対ERSPは同じになり得る。
% ERSPの強さは元EEGへの寄与率ではない。ICの自動採否にも使わない。
% 投影ほぼゼロ/非有限値/不正baselineは0dBで埋めず失敗として記録する。

    config=default_config();
    output_parent_specified=false;
    if nargin>=6 && ~isempty(options)
        if ~isstruct(options) || ~isscalar(options), error('singleIC:Options','optionsはscalar structです。'); end
        names=fieldnames(options);
        output_parent_specified=isfield(options,'output_parent');
        for k=1:numel(names)
            if ~isfield(config,names{k}), error('singleIC:Option','未知の設定: %s',names{k}); end
            config.(names{k})=options.(names{k});
        end
    end
    if nargin>=2 && ~isempty(analysis_channels), config.analysis_channels=analysis_channels; end
    if nargin>=3 && ~isempty(ica_seed), config.ica_seed=ica_seed; end
    if nargin>=4 && ~isempty(ersp_max_db), config.ersp_color_limit_db=ersp_max_db; end
    if nargin<5, ic_numbers=[]; end
    config=validate_config(config,ic_numbers);
    if nargin<1 || isempty(input_root) || strlength(string(input_root))==0
        input_root=uigetdir('','03_手動サッケード、被験者、または被験者一覧を選択');
        if isequal(input_root,0), error('singleIC:Cancelled','中断しました。'); end
    end
    input_root=string(input_root);
    if ~isscalar(input_root) || ismissing(input_root) || ~isfolder(input_root)
        error('singleIC:Input','入力フォルダが存在しません。');
    end
    folders=find_manual_saccade_folders(input_root);
    if isempty(folders), error('singleIC:NoSET','03_手動サッケードが見つかりません。'); end
    if strlength(config.ica_model_file)>0 && numel(folders)~=1
        error('singleIC:ModelScope','既存モデルの指定時は1被験者だけを選択してください。');
    end
    initialize_eeglab(config);
    original_rng=rng; restore_rng=onCleanup(@()rng(original_rng)); %#ok<NASGU>
    program_source=fileread([mfilename('fullpath') '.m']);
    ids=strings(numel(folders),1);
    for s=1:numel(folders), [~,ids(s)]=fileparts(fileparts(folders(s))); end
    if numel(unique(ids))~=numel(ids), error('singleIC:DuplicateSubject','被験者名が重複しています。解析範囲を絞ってください。'); end
    config.output_parent=resolve_output_parent(input_root,folders,config.output_parent,output_parent_specified);
    run_root=new_run_folder(config.output_parent);
    ensure_folder(fullfile(run_root,'images'));
    ensure_folder(fullfile(run_root,'boards'));
    manifest=empty_manifest();
    for s=1:numel(folders)
        subject_folder=string(fileparts(folders(s))); [~,id]=fileparts(subject_folder);
        out=string(fullfile(run_root,'subjects',string(id))); ensure_folder(out);
        subject_manifest=empty_manifest();
        scalp_topography=struct('Values',[],'Locations',[],'StandardLocationLabels',strings(0,1),'Message',"No ICA model");
        fprintf('\n単独IC ERSP [%d/%d] %s\n出力: %s\n',s,numel(folders),id,out);
        try
            datasets=load_and_preprocess_conditions(folders(s),config);
            labels=string({datasets.ConditionLabel});
            if numel(unique(labels))~=numel(labels)
                error('singleIC:DuplicateCondition','同じ条件のSETが複数あります。混ぜずに入力を整理してください。');
            end
            sources=source_inventory(datasets);
            [model,report]=obtain_model(datasets,subject_folder,config,sources);
            scalp_topography=prepare_scalp_topography(model);
            n=size(model.icaweights,1);
            selected=ic_numbers;
            if isempty(selected), selected=1:n; end
            if any(selected>n), error('singleIC:ICRange','この被験者のICは1～%dです。',n); end
            selected=selected(:)';
            icaweights=model.icaweights; icasphere=model.icasphere; icawinv=model.icawinv;
            icachansind=model.icachansind; channel_labels=string({model.chanlocs.labels});
            ica_seed=report.ICASeed; analysis_config=config; %#ok<NASGU>
            model_path=fullfile(out,'共通ICA_単独解析モデル.mat');
            save(model_path,'icaweights','icasphere','icawinv','icachansind','channel_labels', ...
                'ica_seed','analysis_config','report','sources','selected','program_source','scalp_topography','-v7');
            for ic=selected
                ic_folder=fullfile(out,sprintf('IC_%02d',ic)); ensure_folder(ic_folder);
                reference_axes=[];
                for c=1:numel(datasets)
                    d=datasets(c); EEG=project_single_ic(d.EEG,model,ic);
                    try
                        epoch=epoch_for_newtimef(EEG,config);
                        for channel=config.analysis_channels
                            pos=find(strcmpi(config.target_channels,channel),1);
                            ch=d.TargetIndices(pos);
                            row=manifest_row(string(id),ic,d,channel,out,report,config);
                            try
                                weight=model.icawinv(pos,ic);
                                if abs(weight)<=config.minimum_relative_loading*max(abs(model.icawinv(:,ic)))
                                    error('singleIC:ZeroProjection','この電極のIC投影がほぼゼロです。ERSPは解釈できません。');
                                end
                                [trial_ersp,freqs,ersp_mean,times,itc,powbase]= ...
                                    compute_newtimef_trial_ersp(epoch,ch,config);
                                if any(~isfinite(ersp_mean(:))) || ~isreal(ersp_mean)
                                    error('singleIC:InvalidERSP','ERSPに非有限値があります。');
                                end
                                if isempty(reference_axes), reference_axes=struct('Times',times,'Freqs',freqs); end
                                if ~isequal(reference_axes.Times,times) || ~isequal(reference_axes.Freqs,freqs)
                                    error('singleIC:Axes','条件・電極間の時間周波数軸が一致しません。');
                                end
                                [~,ersp_sem,ersp_n]=mean_sem_over_trials(trial_ersp);
                                trial_count=epoch.trials; source_set=d.SourceSET; condition_label=d.ConditionLabel;
                                kept_components=ic; removed_components=setdiff(1:n,ic); scalp_loading=weight;
                                analysis_config=config; ic_report=report; %#ok<NASGU>
                                channel_folder=fullfile(ic_folder,channel); ensure_folder(channel_folder);
                                stem=safe_filename(d.BaseName)+sprintf('_IC%02d_',ic)+channel+'_単独IC_newtimef_ERSP';
                                row.ERSPMAT=string(fullfile(channel_folder,stem+'.mat'));
                                save(row.ERSPMAT,'source_set','condition_label','kept_components','removed_components', ...
                                    'analysis_config','ic_report','scalp_loading','ersp_mean','ersp_sem','ersp_n', ...
                                    'times','freqs','itc','powbase','trial_count','-v7');
                                row.TrialCount=trial_count; row.Status="Complete";
                            catch err
                                row.Status="Failed"; row.Message=string(getReport(err,'extended','hyperlinks','off'));
                                warning('singleIC:CellFailed','IC%d / %s / %s: %s',ic,d.ConditionLabel,channel,err.message);
                            end
                            subject_manifest=[subject_manifest;row]; %#ok<AGROW>
                        end
                    catch err
                        for channel=config.analysis_channels
                            row=manifest_row(string(id),ic,d,channel,out,report,config);
                            row.Status="Failed"; row.Message=string(getReport(err,'extended','hyperlinks','off'));
                            subject_manifest=[subject_manifest;row]; %#ok<AGROW>
                        end
                        warning('singleIC:ConditionFailed','IC%d / %s: %s',ic,d.ConditionLabel,err.message);
                    end
                    writetable(subject_manifest,fullfile(out,'ERSP_manifest.csv'),'Encoding','UTF-8');
                end
                try, render_ic_image(subject_manifest(subject_manifest.IC==ic,:),run_root,string(id),ic,config,scalp_topography);
                catch err, warning('singleIC:Image','%s IC%dの個別画像: %s',id,ic,err.message); end
            end
            save(fullfile(out,'run_metadata.mat'),'config','report','sources','selected','program_source', ...
                'subject_manifest','scalp_topography','-v7');
        catch err
            d=struct('ConditionLabel',"",'SourceSET',folders(s),'RawSaccadeEvents',NaN,'DuplicatesRemoved',NaN);
            row=manifest_row(string(id),NaN,d,"",out,struct('ICASeed',config.ica_seed),config);
            row.Status="Failed"; row.Message=string(getReport(err,'extended','hyperlinks','off'));
            subject_manifest=[subject_manifest;row]; %#ok<AGROW>
            writetable(subject_manifest,fullfile(out,'ERSP_manifest.csv'),'Encoding','UTF-8');
            warning('singleIC:SubjectFailed','%s: %s',id,err.message);
        end
        try, render_subject_board(subject_manifest,run_root,string(id),config,scalp_topography);
        catch err, warning('singleIC:Board','%sの全IC一覧画像: %s',id,err.message); end
        manifest=[manifest;subject_manifest]; %#ok<AGROW>
        writetable(manifest,fullfile(run_root,'ERSP_manifest_all.csv'),'Encoding','UTF-8');
        fprintf('完了セル%d / 失敗セル%d\n',sum(subject_manifest.Status=="Complete"),sum(subject_manifest.Status=="Failed"));
    end
    save(fullfile(run_root,'run_metadata_all.mat'),'config','manifest','program_source','input_root','-v7');
    fprintf('各ICの個別画像: %s\n',fullfile(run_root,'images'));
    fprintf('被験者ごとの全IC一覧画像: %s\n',fullfile(run_root,'boards'));
end
 
function config=default_config()
    config.target_channels=["F3","F4","Fz","O1","O2","Oz","PO7","PO8"];
    config.analysis_channels="Oz";
    config.eog_channel_patterns=["EOG","HEOG","VEOG"];
    config.event_type="Saccade";
    config.event_duplicate_tolerance_samples=1;
    config.analysis_highpass_hz=0.2;
    config.analysis_lowpass_hz=220;
    config.apply_average_reference=true;
    config.ica_highpass_hz=1;
    config.ica_lowpass_hz=100;
    config.ica_training_srate=250;
    config.ica_rank=7;
    config.ica_seed=3;
    config.ersp_color_limit_db=1.5;
    config.display_window_ms=[-200 300];
    config.baseline_window_ms=[-250 -100];
    config.epoch_window_s=[-4 4];
    config.newtimef_frequency_limits_hz=[0.5 200];
    config.ersp_display_max_hz=200;
    config.newtimef_frequency_count=200;
    config.newtimef_frequency_scale='linear';
    config.wavelet_cycles=3;
    config.newtimef_timesout=400;
    config.padratio=1;
    config.ica_model_file="";
    config.output_parent="C:\Users\tatsuya\学校法人東海大学\高雄研　データ共有 - 光瞳孔反射研究プロジェクト\02瞳孔研究2025\05Phantom Array Effect\01　実験データ\06_模擬実験_2610\01_CWT\08_ERSP_QC";
    config.minimum_relative_loading=1e-10;
    config.eeglab_fallback="C:\Users\tatsuya\Downloads\Apps\eeglab2025.1.0";
end

function config=validate_config(config,ics)
    config.target_channels=string(config.target_channels(:)');
    config.analysis_channels=string(config.analysis_channels(:)');
    if isempty(config.target_channels) || numel(unique(lower(config.target_channels)))~=numel(config.target_channels) || ...
            isempty(config.analysis_channels) || numel(unique(lower(config.analysis_channels)))~=numel(config.analysis_channels) || ...
            ~all(ismember(lower(config.analysis_channels),lower(config.target_channels)))
        error('singleIC:Channels','学習/出力電極の設定が不正です。');
    end
    validateattributes(config.ica_seed,{'numeric'},{'scalar','integer','finite','nonnegative','<=',double(intmax('uint32'))});
    validateattributes(config.ica_rank,{'numeric'},{'scalar','integer','finite','positive'});
    validateattributes(config.ersp_color_limit_db,{'numeric'},{'scalar','real','finite','positive'});
    if ~isempty(ics), validateattributes(ics,{'numeric'},{'vector','integer','finite','positive','real'}); end
    if numel(unique(ics))~=numel(ics), error('singleIC:ICDuplicates','IC番号が重複しています。'); end
    for field=["ica_model_file","output_parent","eeglab_fallback"]
        value=string(config.(field));
        if ~isscalar(value) || ismissing(value), error('singleIC:Path','パス設定が不正です: %s',field); end
        config.(field)=value;
    end
    validateattributes(config.apply_average_reference,{'logical'},{'scalar'});
    for field=["epoch_window_s","baseline_window_ms","display_window_ms","newtimef_frequency_limits_hz"]
        validateattributes(config.(field),{'numeric'},{'vector','real','finite','numel',2,'increasing'});
    end
    validateattributes(config.ersp_display_max_hz,{'numeric'}, ...
        {'scalar','real','finite','>',config.newtimef_frequency_limits_hz(1), ...
         '<=',config.newtimef_frequency_limits_hz(2)});
    if config.newtimef_frequency_limits_hz(1)<=0 || config.baseline_window_ms(2)>=0 || ...
            config.baseline_window_ms(1)<1000*config.epoch_window_s(1) || ...
            config.display_window_ms(1)<1000*config.epoch_window_s(1) || config.display_window_ms(2)>1000*config.epoch_window_s(2)
        error('singleIC:Windows','周波数・baseline・表示窓が不正です。');
    end
    for field=["analysis_highpass_hz","analysis_lowpass_hz","ica_highpass_hz","ica_lowpass_hz", ...
            "ica_training_srate","minimum_relative_loading"]
        validateattributes(config.(field),{'numeric'},{'scalar','real','finite','positive'});
    end
    for field=["newtimef_frequency_count","newtimef_timesout","padratio"]
        validateattributes(config.(field),{'numeric'},{'scalar','integer','finite','positive'});
    end
end

function folder=new_run_folder(parent)
    ensure_folder(parent);
    base=string(fullfile(parent,"IC_run_"+string(datetime('now','TimeZone','Asia/Tokyo','Format','yyyyMMdd_HHmmss_SSS'))));
    folder=base; k=0;
    while isfolder(folder), k=k+1; folder=base+"_"+k; end
    ensure_folder(folder);
end

function parent=resolve_output_parent(input_root,manual_folders,configured_parent,explicit_parent)
% 個人選択と「一覧に1名だけある」場合を区別する。明示した出力先は尊重。
    parent=string(configured_parent);
    if explicit_parent && strlength(parent)>0, return; end
    if numel(manual_folders)==1
        manual_folder=string(manual_folders(1));
        subject_folder=string(fileparts(manual_folder));
        selected=string(java.io.File(char(input_root)).getCanonicalPath());
        manual_path=string(java.io.File(char(manual_folder)).getCanonicalPath());
        subject_path=string(java.io.File(char(subject_folder)).getCanonicalPath());
        if strcmpi(selected,manual_path) || strcmpi(selected,subject_path)
            parent=string(fullfile(subject_folder,'08_ERSP_QC'));
            return;
        end
    end
    % 被験者一覧の保存先は従来設定を維持。空指定なら選択範囲の08へ。
    if strlength(parent)==0, parent=string(fullfile(input_root,'08_ERSP_QC')); end
end

function sources=source_inventory(datasets)
    sources=table();
    for k=1:numel(datasets)
        path=datasets(k).SourceSET; file=dir(path); EEG=datasets(k).EEG;
        fdt=dir(replace(path,'.set','.fdt')); fdt_bytes=NaN; fdt_date=NaN;
        if isscalar(fdt), fdt_bytes=fdt.bytes; fdt_date=fdt.datenum; end
        row=table(path,datasets(k).ConditionLabel,double(file.bytes),file.datenum,fdt_bytes,fdt_date,EEG.pnts,EEG.srate, ...
            'VariableNames',{'SourceSET','Condition','SETBytes','SETModifiedDatenum','FDTBytes','FDTModifiedDatenum','Samples','SamplingRate'});
        sources=[sources;row]; %#ok<AGROW>
    end
end

function [model,report]=obtain_model(datasets,subject_folder,config,sources)
    if strlength(config.ica_model_file)==0
        [model,report]=learn_common_model(datasets,config);
        report.ModelOrigin="New common ICA (all conditions, one training)";
        report.ModelSource="";
        return;
    end
    path=config.ica_model_file;
    if ~isfile(path), error('singleIC:ModelFile','ICAモデルがありません: %s',path); end
    saved=load(path);
    required={'icaweights','icasphere','icawinv','icachansind','channel_labels','analysis_config','report'};
    if ~all(isfield(saved,required)), error('singleIC:ModelSchema','共通ICA分解MATを指定してください。'); end
    if isfield(saved,'sources')
        if ~isequaln(saved.sources,sources), error('singleIC:ModelSources','モデル保存時と入力SET/FDT情報が一致しません。'); end
    else
        canonical=char(java.io.File(char(path)).getCanonicalPath());
        subject=char(java.io.File(char(subject_folder)).getCanonicalPath());
        if ~startsWith(lower(string(canonical)),lower(string(subject)+filesep))
            error('singleIC:ModelParticipant','入力照合情報のない旧モデルは、選択被験者フォルダ内のモデルに限定します。');
        end
        warning('singleIC:UnverifiedModelInputs','旧モデルには入力SET照合情報がありません。元解析と同じSETであることを確認してください。');
    end
    for field=["target_channels","analysis_highpass_hz","analysis_lowpass_hz","apply_average_reference", ...
            "ica_highpass_hz","ica_lowpass_hz","ica_training_srate","ica_rank"]
        if ~isfield(saved.analysis_config,field) || ~isequaln(saved.analysis_config.(field),config.(field))
            error('singleIC:ModelPreprocessing','旧モデルと設定が一致しません: %s。optionsで合わせてください。',field);
        end
    end
    labels=string({datasets(1).EEG.chanlocs.labels});
    if ~isequal(lower(string(saved.channel_labels(:)')),lower(labels)) || ...
            ~isequal(double(saved.icachansind(:)'),datasets(1).TargetIndices)
        error('singleIC:ModelChannels','モデルの電極順が入力と一致しません。');
    end
    model=struct('icaweights',saved.icaweights,'icasphere',saved.icasphere, ...
        'icawinv',saved.icawinv,'icachansind',saved.icachansind,'chanlocs',datasets(1).EEG.chanlocs);
    n=size(model.icaweights,1); m=numel(model.icachansind);
    if size(model.icasphere,2)~=m || size(model.icaweights,2)~=size(model.icasphere,1) || ...
            ~isequal(size(model.icawinv),[m n]) || any(~isfinite([model.icaweights(:);model.icasphere(:);model.icawinv(:)]))
        error('singleIC:ModelShape','ICA行列が不正です。');
    end
    if norm(model.icawinv-pinv(model.icaweights*model.icasphere),'fro')>1e-6*max(1,norm(model.icawinv,'fro'))
        error('singleIC:ModelInverse','IC逆投影行列が分解と一致しません。');
    end
    report=saved.report;
    if ~all(isfield(report,{'ICASeed','Classes','Classification'})) || size(report.Classification,1)~=n
        error('singleIC:ModelReport','モデルのICLabel情報が不正です。');
    end
    report.ModelOrigin="Reused saved common ICA (no retraining)"; report.ModelSource=path;
    fprintf('既存ICAを再利用: %s / 保存時seed=%u（今回のseed指定は学習に使わない）\n',path,report.ICASeed);
end

function EEG=project_single_ic(EEG,model,ic)
    indices=model.icachansind;
    activity=double(model.icaweights(ic,:)*model.icasphere)*double(EEG.data(indices,:));
    EEG.data(indices,:)=double(model.icawinv(:,ic))*activity;
    EEG=clear_ica_fields(EEG);
    EEG.etc.single_ic=struct('OriginalIC',ic,'Method',"One-component scalp backprojection, no residual", ...
        'ICATargetIndices',indices,'NonICATargetChannels',"Unchanged preprocessed signals");
    EEG=eeg_checkset(EEG,'eventconsistency');
end

function table_out=empty_manifest()
    table_out=table('Size',[0 14],'VariableTypes',{'string','double','string','string','string', ...
        'double','double','double','double','double','string','string','string','string'}, ...
        'VariableNames',{'ParticipantID','IC','Condition','Channel','SourceSET','ICASeed','ERSPDisplayMaxDb', ...
        'RawSaccadeEvents','DuplicateEventsRemoved','TrialCount','Status','Message','OutputRoot','ERSPMAT'});
end

function row=manifest_row(id,ic,d,channel,out,report,config)
    template=empty_manifest();
    row=table(id,double(ic),string(d.ConditionLabel),string(channel),string(d.SourceSET), ...
        double(report.ICASeed),config.ersp_color_limit_db,double(d.RawSaccadeEvents),double(d.DuplicatesRemoved), ...
        NaN,"Pending","",string(out),"",'VariableNames',template.Properties.VariableNames);
end

function render_ic_image(rows,run_root,id,ic,config,scalp_topography)
    conditions=unique(rows.Condition,'stable'); channels=config.analysis_channels;
    ncols=numel(conditions)+1;
    fig=figure('Visible','off','Color','w','Position',[25 25 max(1150,350+390*numel(conditions)) max(400,240*numel(channels))]);
    cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
    layout=tiledlayout(fig,numel(channels),ncols,'TileSpacing','compact','Padding','compact');
    scalp_ax=nexttile(layout,1,[numel(channels) 1]); render_scalp_map(scalp_ax,scalp_topography,ic);
    for channel_index=1:numel(channels)
        ch=channels(channel_index);
        for condition_index=1:numel(conditions)
            condition=conditions(condition_index);
            ax=nexttile(layout,(channel_index-1)*ncols+condition_index+1);
            selected=rows(rows.Channel==ch & rows.Condition==condition,:);
            if height(selected)==1 && selected.Status=="Complete" && isfile(selected.ERSPMAT)
                data=load(selected.ERSPMAT,'times','freqs','ersp_mean','trial_count');
                imagesc(ax,data.times,data.freqs,data.ersp_mean); set(ax,'YDir','normal');
                xlim(ax,config.display_window_ms);
                ylim(ax,[config.newtimef_frequency_limits_hz(1) config.ersp_display_max_hz]);
                clim(ax,[-config.ersp_color_limit_db config.ersp_color_limit_db]);
                xline(ax,0,'--m'); title(ax,sprintf('%s / %s / N=%d',condition,ch,data.trial_count),'Interpreter','none');
                xlabel(ax,'Time (ms)'); ylabel(ax,'Frequency (Hz)');
            else
                axis(ax,'off'); text(ax,0.5,0.5,"Failed / Missing",'HorizontalAlignment','center','Color',[0.8 0 0]);
                title(ax,condition+" / "+ch,'Interpreter','none');
            end
            colormap(ax,turbo); set(ax,'FontSize',9);
        end
    end
    cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='ERSP (dB)';
    title(layout,sprintf('%s / IC%d ONLY / display %g-%g Hz, +/-%.2f dB', ...
        id,ic,config.newtimef_frequency_limits_hz(1),config.ersp_display_max_hz, ...
        config.ersp_color_limit_db),'Interpreter','none');
    path=fullfile(run_root,'images',safe_filename(id)+sprintf('_IC%02d.png',ic));
    apply_figure_colormaps(fig);
    print_figure_png(fig,path,160);
end

function render_subject_board(rows,run_root,id,config,scalp_topography)
    % MATの同じ解析値を再描画する。個別PNGの縮小・ERSP補間はしない。
    components=unique(rows.IC(isfinite(rows.IC)),'stable');
    conditions=unique(rows.Condition,'stable'); channels=config.analysis_channels;
    if isempty(components)
        fig=figure('Visible','off','Color','w','Position',[25 25 850 400]);
        cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
        ax=axes(fig); axis(ax,'off');
        text(ax,0.5,0.5,'Failed / No IC results','HorizontalAlignment','center','Color',[0.8 0 0]);
        title(ax,id,'Interpreter','none');
    else
        nrows=numel(components)*numel(channels); ncols=numel(conditions)+1;
        % 多電極指定時も省略しない。大きな一覧の詳細はimagesの個別画像で確認。
        fig=figure('Visible','off','Color','w','Position', ...
            [25 25 max(1150,350+390*numel(conditions)) max(400,min(10000,230*nrows))]);
        cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
        layout=tiledlayout(fig,nrows,ncols,'TileSpacing','compact','Padding','compact');
        for ic_index=1:numel(components)
            ic=components(ic_index); first_row=(ic_index-1)*numel(channels);
            scalp_ax=nexttile(layout,first_row*ncols+1,[numel(channels) 1]);
            render_scalp_map(scalp_ax,scalp_topography,ic);
            for channel_index=1:numel(channels)
                ch=channels(channel_index);
                for condition_index=1:numel(conditions)
                    condition=conditions(condition_index);
                    ax=nexttile(layout,(first_row+channel_index-1)*ncols+condition_index+1);
                    selected=rows(rows.IC==ic & rows.Channel==ch & rows.Condition==condition,:);
                    if height(selected)==1 && selected.Status=="Complete" && isfile(selected.ERSPMAT)
                        data=load(selected.ERSPMAT,'times','freqs','ersp_mean','trial_count');
                        imagesc(ax,data.times,data.freqs,data.ersp_mean); set(ax,'YDir','normal');
                        xlim(ax,config.display_window_ms);
                        ylim(ax,[config.newtimef_frequency_limits_hz(1) config.ersp_display_max_hz]);
                        clim(ax,[-config.ersp_color_limit_db config.ersp_color_limit_db]);
                        xline(ax,0,'--m');
                        title(ax,sprintf('IC%d / %s / %s / N=%d',ic,condition,ch,data.trial_count),'Interpreter','none');
                        xlabel(ax,'Time (ms)'); ylabel(ax,'Frequency (Hz)');
                    else
                        axis(ax,'off');
                        text(ax,0.5,0.5,'Failed / Missing','HorizontalAlignment','center','Color',[0.8 0 0]);
                        title(ax,sprintf('IC%d / %s / %s',ic,condition,ch),'Interpreter','none');
                    end
                    colormap(ax,turbo); set(ax,'FontSize',9);
                end
            end
        end
        cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='ERSP (dB)';
        title(layout,sprintf('%s / %d single-IC results / display %g-%g Hz, +/-%.2f dB', ...
            id,numel(components),config.newtimef_frequency_limits_hz(1),config.ersp_display_max_hz, ...
            config.ersp_color_limit_db),'Interpreter','none');
    end
    apply_figure_colormaps(fig);
    print_figure_png(fig,fullfile(run_root,'boards',safe_filename(id)+'.png'),160);
end

function scalp=prepare_scalp_topography(model)
    % 表示専用のコピー。icawinvとicachansindの順番を厳密に対応させる。
    scalp=struct('Values',model.icawinv,'Locations',[], ...
        'StandardLocationLabels',strings(0,1),'Message',"");
    try
        original=model.chanlocs(model.icachansind);
        locations=repmat(struct('labels','','theta',[],'radius',[]),size(original));
        missing=false(size(original));
        for k=1:numel(original)
            converted=original(k);
            if ~valid_polar_location(converted) && ...
                    (valid_coordinate_fields(converted,{'X','Y','Z'}) || ...
                     valid_coordinate_fields(converted,{'sph_theta','sph_phi','sph_radius'}))
                converted=convertlocs(converted,'auto');
            end
            locations(k).labels=original(k).labels;
            if valid_polar_location(converted)
                locations(k).theta=converted.theta; locations(k).radius=converted.radius;
            else
                missing(k)=true;
            end
        end
        if any(missing)
            lookup=resolve_standard_lookup_file();
            standard=readlocs(char(lookup),'defaultelp','BESA');
            positions=find(missing);
            for k=1:numel(positions)
                position=positions(k);
                match=find(strcmpi({standard.labels},locations(position).labels));
                if numel(match)~=1 || ~valid_polar_location(standard(match))
                    error('singleIC:ScalpLocation','頭皮分布用の電極座標がありません: %s',locations(position).labels);
                end
                locations(position).theta=standard(match).theta;
                locations(position).radius=standard(match).radius;
            end
            scalp.StandardLocationLabels=string({locations(missing).labels});
        end
        % 現行EEGLABのreadlocs/topoplotが要求するX/Y/Z等も表示用に生成。
        scalp.Locations=convertlocs(locations,'topo2all');
    catch err
        scalp.Message=string(getReport(err,'extended','hyperlinks','off'));
        warning('singleIC:ScalpLocation','頭皮分布の座標準備に失敗（ERSPは継続）: %s',err.message);
    end
end

function valid=valid_coordinate_fields(location,names)
    valid=all(isfield(location,names));
    if ~valid, return; end
    for k=1:numel(names)
        value=location.(names{k});
        if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
            valid=false; return;
        end
    end
end

function valid=valid_polar_location(location)
    valid=isfield(location,'theta') && isfield(location,'radius') && ...
        isnumeric(location.theta) && isscalar(location.theta) && isfinite(location.theta) && ...
        isnumeric(location.radius) && isscalar(location.radius) && isfinite(location.radius) && location.radius>=0;
end

function render_scalp_map(ax,scalp,ic)
    set(ax,'Tag','singleICScalp');
    try
        if strlength(scalp.Message)>0 || isempty(scalp.Locations)
            error('singleIC:ScalpUnavailable','%s',scalp.Message);
        end
        weights=scalp.Values(:,ic); scale=max(abs(weights));
        if ~isfinite(scale) || scale<=0, error('singleIC:ScalpZero','IC投影係数が不正です。'); end
        % topoplotはaxesのPositionとfigure全体のcolormapを変更するため、
        % 独立した非表示図で描画し、描画オブジェクトだけを一覧のaxesへ移す。
        % ERSP値/頭皮分布を画像化して補間・縮小する処理ではない。
        temp_fig=figure('Visible','off','Color','w','Position',[50 50 450 450]);
        close_temp=onCleanup(@()close(temp_fig)); %#ok<NASGU>
        temp_ax=axes(temp_fig);
        previous_warning=warning;
        restore_warning=onCleanup(@()warning(previous_warning)); %#ok<NASGU>
        topoplot(weights/scale,scalp.Locations,'maplimits',[-1 1], ...
            'electrodes','labels','whitebk','on','verbose','off');
        clear restore_warning;
        copyobj(get(temp_ax,'Children'),ax);
        set(ax,'XLim',get(temp_ax,'XLim'),'YLim',get(temp_ax,'YLim'), ...
            'ZLim',get(temp_ax,'ZLim'),'View',get(temp_ax,'View'));
        axis(ax,'equal'); axis(ax,'off');
        colormap(ax,jet(256)); clim(ax,[-1 1]);
        cb=colorbar(ax); cb.Label.String='Normalized loading (not dB)';
        cb.Ticks=[-1 0 1]; set(ax,'FontSize',9);
        heading=sprintf('IC%d / Scalp map',ic);
        if ~isempty(scalp.StandardLocationLabels), heading=[heading ' (standard locs)']; end
        title(ax,heading,'Interpreter','none');
    catch err
        cla(ax); axis(ax,'off');
        title(ax,sprintf('IC%d / Scalp map',ic),'Interpreter','none');
        text(ax,0.5,0.5,'Scalp map unavailable','HorizontalAlignment','center','Color',[0.8 0 0]);
        warning('singleIC:ScalpPlot','IC%dの頭皮分布を表示できません（ERSPは継続）: %s',ic,err.message);
    end
end

function apply_figure_colormaps(fig)
    % topoplot内部のfigure共通colormap変更後に、全axesを明示して固定する。
    for ax=findall(fig,'Type','axes')'
        if strcmp(get(ax,'Tag'),'singleICScalp'), colormap(ax,jet(256));
        else, colormap(ax,turbo); end
    end
end

% 以下の前処理・連結・newtimefヘルパーはCWT_20261006の2026-10-07版から移植。
function [model, report] = learn_common_model(datasets, config)
    [merged, ~, ~] = ...
        concatenate_condition_data(datasets);
    target_indices = datasets(1).TargetIndices;

    try
        if ~isfield(merged.chanlocs, 'X') || ...
                isempty(merged.chanlocs(target_indices(1)).X)
            lookup_file = resolve_standard_lookup_file();
            merged = pop_chanedit(merged, ...
                'lookup', char(lookup_file));
        end
    catch ME
        error('singleIC:ChannelLookup', ...
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
    fprintf('  boundary保持: 連結時%d個 / ICA学習時%d個（条件接続点%d個）\n', ...
        numel(merged.event), numel(training.event), numel(datasets) - 1);

    rng(config.ica_seed, 'twister');
    data_rank = rank(double(training.data(target_indices, :)));
    rank_value = min([config.ica_rank, data_rank, ...
        numel(target_indices)]);
    if rank_value<1, error('singleIC:Rank','ICA学習データのrankが0です。'); end
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

    classes=["Brain","Muscle","Eye","Heart","Line Noise","Channel Noise","Other"];
    classification=nan(size(merged.icaweights,1),numel(classes));
    report=struct('Rank',rank_value,'DataRank',data_rank,'RequestedRank',config.ica_rank, ...
        'ICASeed',config.ica_seed,'Classes',classes,'Classification',classification, ...
        'ICATrainingSrate',training.srate,'AnalysisFilterHz',[config.analysis_highpass_hz config.analysis_lowpass_hz], ...
        'BoundaryHandling',"Preserve boundaries; hard boundaries between conditions", ...
        'AnalysisBoundaryCount',numel(merged.event),'ICATrainingBoundaryCount',numel(training.event), ...
        'SelectionMethod',"Each IC separately; no ICLabel selection");
    model=struct('icaweights',merged.icaweights,'icasphere',merged.icasphere, ...
        'icawinv',merged.icawinv,'icachansind',merged.icachansind,'chanlocs',merged.chanlocs);
end

function initialize_eeglab(config)
    if exist('eeglab', 'file') ~= 2 && isfolder(config.eeglab_fallback)
        addpath(config.eeglab_fallback);
    end
    if exist('eeglab', 'file') ~= 2
        error('singleIC:EEGLABNotFound', ...
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
    if isempty(set_list)
        error('singleIC:TooFewConditions', ...
            'SETがありません: %s', manual_folder);
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
            error('singleIC:NotContinuous', ...
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
            error('singleIC:ChannelMismatch', ...
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
        % 元のSET名・SourceSETは変えず、生成物のベース名だけ正規化。
        output_base_name = string(base_name);
        if is_green
            output_base_name = replace(output_base_name, "緑なし", "NoGreen");
        end
        datasets(file_index).BaseName = output_base_name;
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

function ensure_folder(folder_path)
    if ~isfolder(folder_path)
        [created, message] = mkdir(folder_path);
        if ~created
            error('singleIC:CreateFolderFailed', ...
                '出力フォルダを作成できません: %s (%s)', ...
                folder_path, message);
        end
    end
end

function lookup_file = resolve_standard_lookup_file()
    eeglab_file = string(which('eeglab'));
    if strlength(eeglab_file) == 0
        error('singleIC:EEGLABPathMissing', ...
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
        error('singleIC:LookupFileMissing', ...
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
                error('singleIC:FDTNotFound', ...
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
    error('singleIC:NonfiniteData', ...
        ['非有限値を含む%dサンプルがあります。線形補間は行わない設定のため、' ...
         '入力データを確認してください。'], sum(invalid));
end

function target_indices = find_target_channels(EEG, target_channels)
    labels = string({EEG.chanlocs.labels});
    target_indices = zeros(1, numel(target_channels));
    for channel_index = 1:numel(target_channels)
        found = find(strcmpi(labels, target_channels(channel_index)), 1);
        if isempty(found)
            error('singleIC:MissingChannel', ...
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
    error('singleIC:EOGNotFound', ...
        'EOGチャンネルが見つかりません。候補: %s', ...
        strjoin(patterns, ', '));
end

function [label, value, is_green] = parse_condition_name(base_name)
    base_name = string(base_name);
    is_green = contains(base_name, "緑なし") || ...
        contains(base_name, "NoGreen", 'IgnoreCase', true);
    if is_green
        label = "NoGreen";
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
            error('singleIC:DatasetMismatch', ...
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
    % 連結ICA用には不連続点だけを引き継ぐ。Saccade等の解析イベントは
    % 元の各datasetsに残し、ICA適用後の条件別エポック化でそのまま使う。
    merged.event = struct('type', {}, 'latency', {}, 'duration', {});
    for data_index = 1:n_datasets
        EEG = datasets(data_index).EEG;
        types = event_types_as_strings(EEG.event);
        boundary_indices = union(find(strcmpi(types, "boundary")), ...
            eeg_findboundaries(EEG));
        for event_index = boundary_indices(:)'
            source_event = EEG.event(event_index);
            latency = double(source_event.latency);
            if ~isscalar(latency) || ~isfinite(latency) || ...
                    latency < 0.5 || latency > EEG.pnts + 0.5
                error('singleIC:InvalidBoundaryLatency', ...
                    '条件%dのboundary位置が不正です。入力SETを確認してください。', ...
                    data_index);
            end
            duration = 0;
            if isfield(source_event, 'duration') && ...
                    ~isempty(source_event.duration)
                duration = double(source_event.duration);
            end
            if ~isscalar(duration) || ...
                    ~(isnan(duration) || (isfinite(duration) && duration >= 0))
                error('singleIC:InvalidBoundaryDuration', ...
                    '条件%dのboundary durationが不正です。', data_index);
            end
            merged.event(end + 1) = struct('type', 'boundary', ...
                'latency', latency + segment_starts(data_index) - 1, ...
                'duration', duration);
        end
        if data_index > 1
            % EEGLAB pop_mergesetと同じ半サンプル位置・NaN duration。
            % 実際のデータを追加・削除せず、接続点を不連続として区切る。
            merged.event(end + 1) = struct('type', 'boundary', ...
                'latency', segment_starts(data_index) - 0.5, ...
                'duration', NaN);
        end
    end
    if ~isempty(merged.event)
        [~, event_order] = sort([merged.event.latency]);
        merged.event = merged.event(event_order);
    end
    merged.urevent = struct([]);
    merged.epoch = [];
    merged.setname = '全条件連結_共通ICA';
    merged = clear_ica_fields(merged);
    merged = eeg_checkset(merged, 'eventconsistency');
end

function EEG_epoch = epoch_for_newtimef(EEG, config)
    EEG_epoch = pop_epoch(EEG, {char(config.event_type)}, ...
        config.epoch_window_s, 'epochinfo', 'yes');
    if EEG_epoch.trials < 2
        error('singleIC:TooFewEpochs', ...
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
        error('singleIC:InvalidBaseline', ...
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
    error('singleIC:UnexpectedTFSize', ...
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

function safe_name = safe_filename(input_name)
    safe_name = regexprep(string(input_name), '[<>:"/\\|?*]', '_');
    safe_name = strip(safe_name);
    if strlength(safe_name) == 0
        safe_name = "unnamed";
    end
end

function print_figure_png(fig, output_png, resolution)
    % exportgraphicsを多数回呼び出した際のgraphics handshaking timeoutを
    % 避けるため、非表示figureは従来型のprintでPNGへ保存する。
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    print(fig, char(output_png), '-dpng', ...
        sprintf('-r%d', resolution));
end
