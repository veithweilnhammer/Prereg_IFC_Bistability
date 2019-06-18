%% Coordination for batch SPM analysis: PE_bistability%%
clear all
close all

spm('defaults', 'FMRI');
rootdir = '/media/veithweilnhammer/My Passport/PE_Bistability/';
cd(rootdir)

create_folder_structure = 0;
import_dicoms = 0;
preprocess_data = 0;
prepare_behavioural_data = 0;
estimate_GLM = 0;
build_contrasts = 0;
correct_contrasts = 0;
define_ANOVA = 1;
estimate_ANOVA = 1;

subjects = {'PE_observer_001'; 'PE_observer_005'; 'PE_observer_010';};
%subjects = {'PE_observer_011'};

run_type = 'run' %%;  run

experiments_to_preprocess = {'control'}


for i = 1:length(subjects) %% loop over participants
    
    parent_folder = fullfile(rootdir,subjects{i});
    %% create folders for analysis
    if create_folder_structure
        
        subfolders = {'anatomy'; 'anatomy/localizer/'; 'dicoms'; 'fieldmaps'; 'FLAIR';...
            'functional'; 'functional/ambiguity';'functional/graded_ambiguity'; 'functional/control'; ...
            'response'; 'response/ambiguity';'response/graded_ambiguity'; 'response/control'; 'response/modeling'; ...
            'response/design_files'; 'stats'; 'stats/ambiguity'; 'stats/graded_ambiguity';...
            'stats/control';'batches'; 'report_quality'};
        
        
        for ii = 1:length(subfolders)
            if ~exist(fullfile(parent_folder, subfolders{ii}))
                mkdir(fullfile(parent_folder, subfolders{ii}));
            end
        end
    end
    
    %% import dicoms into approrpiate folders
    
    if import_dicoms
        
        folders_to_convert = dir(fullfile(parent_folder,'dicoms','0*'))
        
        %% Check for every participant!
        
        % Participant
        Localizer{1} = '001_AAHead_Scout_64ch-head-coil';
        Localizer{2} = '002_AAHead_Scout_64ch-head-coil';
        Localizer{3} = '003_AAHead_Scout_64ch-head-coil';
        Localizer{4} = '004_AAHead_Scout_64ch-head-coil';
        
        Fieldmap{1} = '007_gre_field_mapping';
        Fieldmap{2} = '008_gre_field_mapping';
        
        FLAIR = '006_t2_tse_dark-fluid_fs_tra';
        
        Anatomy='005_MPRAGE_freesurfer';
        
        Ambiguity_Run{1} = '009_ep2d_bold';
        
        Graded_Ambiguity_Run{1} = '010_ep2d_bold';
        Graded_Ambiguity_Run{2} = '011_ep2d_bold';
        Graded_Ambiguity_Run{3} = '012_ep2d_bold';
        
        Control_Run{1} = '013_ep2d_bold';
        
        clear matlabbatch
        spm_jobman('initcfg');
        
        job_n=0;
        for ii = 1:length(folders_to_convert)
            
            outdir = [];
            
            if strcmp(Localizer{1}, folders_to_convert(ii).name) == 1 | strcmp(Localizer{2}, folders_to_convert(ii).name) == 1 ...
                    | strcmp(Localizer{3}, folders_to_convert(ii).name) == 1 | strcmp(Localizer{4}, folders_to_convert(ii).name) == 1
                outdir = fullfile(parent_folder, 'anatomy','localizer')
                
            elseif strcmp(Fieldmap{1}, folders_to_convert(ii).name) == 1
                outdir = fullfile(parent_folder, 'fieldmaps','fieldmapAP')
            elseif strcmp(Fieldmap{2}, folders_to_convert(ii).name) == 1
                outdir = fullfile(parent_folder, 'fieldmaps','fieldmapPA')
                
            elseif strcmp(FLAIR, folders_to_convert(ii).name) == 1
                outdir = fullfile(parent_folder, 'FLAIR')
                
            elseif strcmp(Anatomy, folders_to_convert(ii).name) == 1
                outdir = fullfile(parent_folder, 'anatomy')
                
            else
                for idx = 1 : length(Ambiguity_Run)
                    if strcmp(Ambiguity_Run{idx}, folders_to_convert(ii).name) == 1
                        outdir = fullfile(parent_folder, 'functional','ambiguity',['run' num2str(idx)])
                    end
                end
                
                for idx = 1 : length(Graded_Ambiguity_Run)
                    if strcmp(Graded_Ambiguity_Run{idx}, folders_to_convert(ii).name) == 1
                        outdir = fullfile(parent_folder, 'functional','graded_ambiguity',['run' num2str(idx)])
                    end
                end
                
                for idx = 1 : length(Control_Run)
                    if strcmp(Control_Run{idx}, folders_to_convert(ii).name) == 1
                        outdir = fullfile(parent_folder, 'functional','control',['run' num2str(idx)])
                    end
                end
                
            end
            
            if ~isempty(outdir)
                job_n = job_n + 1;
                
                files=dir(fullfile(parent_folder,'dicoms',folders_to_convert(ii).name, '*dcm'));
                for f=1:length(files)
                    matlabbatch{job_n}.spm.util.import.dicom.data(f,1) = cellstr(fullfile(parent_folder,'dicoms',folders_to_convert(ii).name, files(f).name));
                end
                matlabbatch{job_n}.spm.util.import.dicom.root = 'flat';
                
                if ~exist(outdir)
                    mkdir(outdir);
                end
                
                matlabbatch{job_n}.spm.util.import.dicom.outdir = cellstr(outdir);
                matlabbatch{job_n}.spm.util.import.dicom.protfilter = '.*';
                matlabbatch{job_n}.spm.util.import.dicom.convopts.format = 'nii';
                matlabbatch{job_n}.spm.util.import.dicom.convopts.icedims = 0;
                
            end
            
        end
        
        spm_jobman(run_type, matlabbatch);
        spm('show');
    end
    
    %% preprocess data
    
    
    if preprocess_data
        clear matlabbatch
        
        for exp = 1:length(experiments_to_preprocess)
            
            run_folders = dir(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, 'run*'))
            
            preproc_branches = [8];
            
            clear matlabbatch
            spm_jobman('initcfg');
            
            for ii = 1:length(run_folders)
                if ~exist(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, run_folders(ii).name, [num2str(preproc_branches) 'mm']))
                    mkdir(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, run_folders(ii).name, [num2str(preproc_branches) 'mm']));
                end
            end
            
            cd(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}))
            preproc_fmri;
            save(fullfile(parent_folder, 'batches', [experiments_to_preprocess{exp} 'preproc_batch.mat']), 'matlabbatch')
            
            
            spm_jobman(run_type, matlabbatch)
            
            %% move preprocessed files to appropriate folder
            for ii=1:length(run_folders)
                files_to_move = dir(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, run_folders(ii).name, '*af*'));
                for move_it = 1:length(files_to_move)
                    movefile(fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, run_folders(ii).name, files_to_move(move_it).name), ...
                        fullfile(parent_folder, 'functional', experiments_to_preprocess{exp}, run_folders(ii).name, [num2str(preproc_branches) 'mm'], files_to_move(move_it).name))
                end
            end
        end
    end
    
    
    if prepare_behavioural_data
        
        if strcmp(experiments_to_preprocess, 'ambiguity') == 1
            clear matlabbatch
            response_files = dir(fullfile(parent_folder, 'response', 'ambiguity', 'Results*'))
            
            %% make folder for stats
            folder_for_design = fullfile(parent_folder, 'stats', 'ambiguity', 'Conventional');
            
            if ~exist(folder_for_design)
                mkdir(folder_for_design);
            end
            
            
            %% general settings for SPM
            matlabbatch{1}.spm.stats.fmri_spec.dir = {folder_for_design};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            
            
            
            for idx = 1:length(response_files)
                
                clear Results
                load(fullfile(parent_folder, 'response', 'ambiguity', response_files(idx).name))
                clear T B
                for trial = 1:length(Results.PDir)
                    
                    T{trial} = [];
                    T{trial}(1,:) = Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                    T{trial}(2,:) = zeros(length(Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime),1);
                    
                    if size(T{trial},2) <= 1
                        B{trial}(1,:) = Results.discrete{trial}(1).*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                    else
                        B{trial}= [];
                        
                    end
                end
                
                %% get motion parameters
                motion_parameters = dir(fullfile(parent_folder, 'functional', 'ambiguity', ['run' num2str(idx)], '8mm', '*.txt'));
                
                %% get scans
                clear files
                files=dir(fullfile(parent_folder, 'functional', 'ambiguity', ['run' num2str(idx)], '8mm', 'sw*'));
                
                for f=1:length(files)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).scans(f,1) = cellstr(fullfile(parent_folder, 'functional', 'ambiguity', ['run' num2str(idx)], '8mm', files(f).name));
                end
                
                %% Matlabbatch
                
                %% Transitions
                Transitions = [T{:}];
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(1).name = 'Transition';
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(1).onset = Transitions(1,:)';
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(1).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(1).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(1).orth = 1;
                
                %% No event blocks
                if ~isempty([B{:}])
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(2).name = 'No_Events';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(2).onset = [B{:}]';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(2).duration = 120;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(2).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(2).orth = 1;
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi_reg = {fullfile(parent_folder, 'functional', 'ambiguity', ['run' num2str(idx)], '8mm', motion_parameters.name)};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).hpf = 128;
                
                %% general settings for HRF
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                
            end
            
            spm_jobman(run_type,matlabbatch);
            
        elseif strcmp(experiments_to_preprocess, 'graded_ambiguity') == 1
            clear matlabbatch
            response_files = dir(fullfile(parent_folder, 'response', 'graded_ambiguity', 'Results*'))
            
            %% make folder for stats
            folder_for_design = fullfile(parent_folder, 'stats', 'graded_ambiguity', 'Conventional');
            
            if ~exist(folder_for_design)
                mkdir(folder_for_design);
            end
            
            
            %% general settings for SPM
            matlabbatch{1}.spm.stats.fmri_spec.dir = {folder_for_design};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            
            for idx = 1:length(response_files)
                
                clear Results
                load(fullfile(parent_folder, 'response', 'graded_ambiguity', response_files(idx).name))
                overlap_timing = [1:Results.Stimulus.frames_per_rot/Results.n_overlaps:Results.Stimulus.frames_per_rot*Results.rot_per_trial].*Results.Monitor.ifi-Results.Monitor.ifi;
                clear T B
                
                %% prepare behavioural data
                for trial = 1:length(Results.PDir)
                    
                    T{trial} = [];
                    T{trial}(1,:) = Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                    T{trial}(2,:) = zeros(length(Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime),1);
                    
                    if size(T{trial},2) <= 1
                        B{trial}(1,:) = Results.discrete{trial}(1).*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                        Cong{trial} = [];
                        Incong{trial} = [];
                    else
                        B{trial}= [];
                        
                        clear congruency congruency_change end_of_events
                        congruency = Results.discrete_steps{trial}-Results.template.discrete_steps{trial};
                        congruency_change = find([1 diff(congruency)]~=0);
                        end_of_events = [[overlap_timing(congruency_change) Results.Stimulus.frames_per_rot*Results.rot_per_trial*Results.Monitor.ifi]...
                            + Results.TrialStartTime{trial}-Results.SessionStartTime];
                        
                        %% Congruency
                        Cong{trial} = overlap_timing(congruency_change(congruency(congruency_change) == 0)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset congruency
                        
                        if ~isempty(Cong{trial})
                            % duration congruency
                            for events = 1:size(Cong{trial},2)
                                Cong{trial}(2, events) = min(end_of_events(end_of_events > Cong{trial}(1, events))) - Cong{trial}(1, events);
                            end
                        end
                        
                        %% Incongruency
                        Incong{trial} = overlap_timing(congruency_change(congruency(congruency_change) ~= 0)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset incongruency
                        
                        if ~isempty(Incong{trial})
                            % duration incongruency
                            for events = 1:size(Incong{trial},2)
                                Incong{trial}(2, events) = min(end_of_events(end_of_events > Incong{trial}(1, events))) - Incong{trial}(1, events);
                            end
                        end
                        
                    end
                end
                
                %% get motion parameters
                motion_parameters = dir(fullfile(parent_folder, 'functional', 'graded_ambiguity', ['run' num2str(idx)], '8mm', '*.txt'));
                
                %% get scans
                clear files
                files=dir(fullfile(parent_folder, 'functional', 'graded_ambiguity', ['run' num2str(idx)], '8mm', 'sw*'));
                
                for f=1:length(files)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).scans(f,1) = cellstr(fullfile(parent_folder, 'functional', 'graded_ambiguity', ['run' num2str(idx)], '8mm', files(f).name));
                end
                
                
                names_of_regressors = {'C1D1'; 'C1D2'; 'C1D3'; 'C1D4'; 'C1D5'; 'C1D6'; ...
                    'C2D1'; 'C2D2'; 'C2D3'; 'C2D4'; 'C2D5'; 'C2D6'; ...
                    'T'; 'B'};
                Transitions = [T{:}];
                
                %% prepare matlabbatch for congruency and incongruency
                deleter = [];
                for trial =  1:length(Results.PDir)
                    clear order_of_block Regress Regress_2
                    order_of_block = find(sort(Results.disambiguation) == Results.disambiguation(trial));
                    
                    
                    Regress = Cong{trial};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).name = names_of_regressors{order_of_block};
                    
                    if ~isempty(Regress)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).onset = Regress(1,:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).duration = Regress(2,:);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).onset = [];
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).duration = [];
                        deleter = [deleter order_of_block];
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).orth = 1;
                    
                    
                    Regress_2 = Incong{trial};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).name = names_of_regressors{order_of_block+6};
                    
                    if ~isempty(Regress_2)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).onset = Regress_2(1,:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).duration = Regress_2(2,:);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).onset = [];
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).duration = [];
                        deleter = [deleter order_of_block+6];
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).orth = 1;
                    
                end
                
                %% delete empty entries
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(deleter) = [];
                final_regress = length(matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond);
                
                %% Transitions
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).name = 'Transition';
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).onset = Transitions(1,:);
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).orth = 1;
                
                %% No event blocks
                if ~isempty([B{:}])
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).name = 'No_Events';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).onset = [B{:}]';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).duration = 120;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).orth = 1;
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi_reg = {fullfile(parent_folder, 'functional', 'graded_ambiguity', ['run' num2str(idx)], '8mm', motion_parameters.name)};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).hpf = 128;
                
                %% general settings for HRF
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                
                
                
            end
            spm_jobman(run_type,matlabbatch);
            
            
        elseif strcmp(experiments_to_preprocess, 'control') == 1
            clear matlabbatch
            response_files = dir(fullfile(parent_folder, 'response', 'control', 'Results*'))
            
            %% make folder for stats
            folder_for_design = fullfile(parent_folder, 'stats', 'control', 'Conventional');
            
            if ~exist(folder_for_design)
                mkdir(folder_for_design);
            end
            
            
            %% general settings for SPM
            matlabbatch{1}.spm.stats.fmri_spec.dir = {folder_for_design};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            
            for idx = 1:length(response_files)
                
                clear Results
                load(fullfile(parent_folder, 'response', 'control', response_files(idx).name))
                overlap_timing = [1:Results.Stimulus.frames_per_cycle/Results.n_overlaps:Results.Stimulus.frames_per_cycle*Results.rot_per_trial].*Results.Monitor.ifi-Results.Monitor.ifi;
                clear T B
                
                %% prepare behavioural data
                for trial = 1:length(Results.PDir)
                    
                    T{trial} = [];
                    T{trial}(1,:) = Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                    T{trial}(2,:) = zeros(length(Results.discrete{trial}.*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime),1);
                    
                    if size(T{trial},2) <= 1
                        B{trial}(1,:) = Results.discrete{trial}(1).*Results.Monitor.ifi+Results.TrialStartTime{trial}-Results.SessionStartTime;
                        Longer{trial} = [];
                        Shorter{trial} = [];
                    else
                        B{trial}= [];
                        
                        clear congruency congruency_change end_of_events
                        longer = zeros(1,length(Results.discrete_steps{trial}));
                        longer(find(Results.discrete_steps{trial} == mode(Results.discrete_steps{trial}))) = 1;
                        longer_change = find([1 diff(longer)]~=0);
                        
                        clear congruency congruency_change end_of_events
                        
                        end_of_events = [[overlap_timing(longer_change) Results.Stimulus.frames_per_cycle*Results.rot_per_trial*Results.Monitor.ifi]...
                            + Results.TrialStartTime{trial}-Results.SessionStartTime];
                        
                        %% Congruency
                        Longer{trial} = overlap_timing(longer_change(longer(longer_change) == 1)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset congruency
                        
                        if ~isempty(Longer{trial})
                            % duration congruency
                            for events = 1:size(Longer{trial},2)
                                Longer{trial}(2, events) = min(end_of_events(end_of_events > Longer{trial}(1, events))) - Longer{trial}(1, events);
                            end
                        end
                        
                        %% Incongruency
                        Shorter{trial} = overlap_timing(longer_change(longer(longer_change) == 0)) + Results.TrialStartTime{trial}-Results.SessionStartTime;  % onset incongruency
                        
                        if ~isempty(Shorter{trial})
                            % duration incongruency
                            for events = 1:size(Shorter{trial},2)
                                Shorter{trial}(2, events) = min(end_of_events(end_of_events > Shorter{trial}(1, events))) - Shorter{trial}(1, events);
                            end
                        end
                        
                    end
                end
                
                %% get motion parameters
                motion_parameters = dir(fullfile(parent_folder, 'functional', 'control', ['run' num2str(idx)], '8mm', '*.txt'));
                
                %% get scans
                clear files
                files=dir(fullfile(parent_folder, 'functional', 'control', ['run' num2str(idx)], '8mm', 'sw*'));
                
                for f=1:length(files)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).scans(f,1) = cellstr(fullfile(parent_folder, 'functional', 'control', ['run' num2str(idx)], '8mm', files(f).name));
                end
                
                
                names_of_regressors = {'A1I1'; 'A1I2'; 'A1I3'; 'A1I4'; 'A1I5'; 'A1I6'; ...
                    'A2I1'; 'A2I2'; 'A2I3'; 'A2I4'; 'A2I5'; 'A2I6'; ...
                    'T'; 'B'};
                Transitions = [T{:}];
                
                %% prepare matlabbatch for longer and shorter
                deleter = [];
                for trial =  1:length(Results.PDir)
                    clear order_of_block Regress Regress_2
                    if length(unique(Results.balance)) < length(Results.balance)
                    end
                    order_of_block = find(sort(Results.balance) == Results.balance(trial));
                    
                    
                    Regress = Longer{trial};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).name = names_of_regressors{order_of_block};
                    
                    if ~isempty(Regress)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).onset = Regress(1,:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).duration = Regress(2,:);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).onset = [];
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).duration = [];
                        deleter = [deleter order_of_block];
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block).orth = 1;
                    
                    
                    Regress_2 = Shorter{trial};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).name = names_of_regressors{order_of_block+6};
                    
                    if ~isempty(Regress_2)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).onset = Regress_2(1,:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).duration = Regress_2(2,:);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).onset = [];
                        matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).duration = [];
                        deleter = [deleter order_of_block+6];
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(order_of_block+6).orth = 1;
                    
                end
                
                %% delete empty entries
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(deleter) = [];
                final_regress = length(matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond);
                
                %% Transitions
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).name = 'Transition';
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).onset = Transitions(1,:);
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+1).orth = 1;
                
                %% No event blocks
                if ~isempty([B{:}])
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).name = 'No_Events';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).onset = [B{:}]';
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).duration = 120;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(idx).cond(final_regress+2).orth = 1;
                end
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).regress = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).multi_reg = {fullfile(parent_folder, 'functional', 'control', ['run' num2str(idx)], '8mm', motion_parameters.name)};
                matlabbatch{1}.spm.stats.fmri_spec.sess(idx).hpf = 128;
                
                %% general settings for HRF
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                
                
                
            end
            spm_jobman(run_type,matlabbatch);
        end
        
    end
    
    if estimate_GLM
        folder_for_estimation = fullfile(parent_folder, 'stats', experiments_to_preprocess, 'Conventional'); % = fullfile(parent_folder, 'stats', 'main_experiment', which_preproc, 'individual_model_based');
        
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile(folder_for_estimation, 'SPM.mat');
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        spm_jobman(run_type,matlabbatch);
    end
    
    if build_contrasts
        
        if strcmp(experiments_to_preprocess, 'ambiguity') == 1
            
        elseif strcmp(experiments_to_preprocess, 'graded_ambiguity') == 1
            
            clear matlabbatch SPM
            
            folder_for_contrast = fullfile(parent_folder, 'stats', 'graded_ambiguity', 'Conventional');
            
            load(fullfile(folder_for_contrast,'SPM.mat'));
            all_betas = SPM.xX.name;
            
            for idx = 1 : length(all_betas)
                all_betas{idx} = all_betas{idx}(7:end);
            end
            
            names_of_regressors = {'C1D1*bf(1)'; 'C1D2*bf(1)'; 'C1D3*bf(1)'; 'C1D4*bf(1)'; 'C1D5*bf(1)'; 'C1D6*bf(1)'; ...
                'C2D1*bf(1)'; 'C2D2*bf(1)'; 'C2D3*bf(1)'; 'C2D4*bf(1)'; 'C2D5*bf(1)'; 'C2D6*bf(1)'; ...
                'Transition*bf(1)'; 'No_Events*bf(1)'};
            
            
            %% make contrast vectors
            n_contrast = 0; clear C
            for idx = 1:length(names_of_regressors)
                for pos_contrast = 1:length(all_betas)
                    
                    if strfind(names_of_regressors{idx}, all_betas{pos_contrast})
                        C(idx,pos_contrast) =  1;
                    else
                        C(idx,pos_contrast) =  0;
                    end
                    
                    
                end
                if any(C(idx,:))
                    n_contrast = n_contrast +1;
                    
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = names_of_regressors{idx};
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = C(idx,:);
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
                end
            end
            
            n_contrast = n_contrast + 1;
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = 'Incongruency_vs_Congruency';
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = sum(C(7:12,:))-sum(C(1:6,:));
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
            
            n_contrast = n_contrast + 1;
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = 'Congruency_vs_Incongruency';
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = -sum(C(7:12,:))+sum(C(1:6,:));
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
            
            
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(folder_for_contrast,'SPM.mat') };
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            spm_jobman(run_type,matlabbatch);
            
            
        elseif strcmp(experiments_to_preprocess, 'control') == 1
            
            clear matlabbatch SPM
            
            folder_for_contrast = fullfile(parent_folder, 'stats', 'control', 'Conventional');
            
            load(fullfile(folder_for_contrast,'SPM.mat'));
            all_betas = SPM.xX.name;
            
            for idx = 1 : length(all_betas)
                all_betas{idx} = all_betas{idx}(7:end);
            end
            
            names_of_regressors = {'A1I1*bf(1)'; 'A1I2*bf(1)'; 'A1I3*bf(1)'; 'A1I4*bf(1)'; 'A1I5*bf(1)'; 'A1I6*bf(1)'; ...
                'A2I1*bf(1)'; 'A2I2*bf(1)'; 'A2I3*bf(1)'; 'A2I4*bf(1)'; 'A2I5*bf(1)'; 'A2I6*bf(1)'; ...
                'Transition*bf(1)'; 'No_Events*bf(1)'};
            
            
            %% make contrast vectors
            n_contrast = 0; clear C
            for idx = 1:length(names_of_regressors)
                for pos_contrast = 1:length(all_betas)
                    
                    if strfind(names_of_regressors{idx}, all_betas{pos_contrast})
                        C(idx,pos_contrast) =  1;
                    else
                        C(idx,pos_contrast) =  0;
                    end
                    
                    
                end
                if any(C(idx,:))
                    n_contrast = n_contrast +1;
                    
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = names_of_regressors{idx};
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = C(idx,:);
                    matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
                end
            end
            
            n_contrast = n_contrast + 1;
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = 'Shorter_vs_Longer';
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = sum(C(7:12,:))-sum(C(1:6,:));
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
            
            n_contrast = n_contrast + 1;
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.name = 'Longer_vs_Shorter';
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.weights = -sum(C(7:12,:))+sum(C(1:6,:));
            matlabbatch{1}.spm.stats.con.consess{n_contrast}.tcon.sessrep = 'none';
            
            
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(folder_for_contrast,'SPM.mat') };
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            spm_jobman(run_type,matlabbatch);
            
            
        end
        
    end
    
    
    if correct_contrasts
        
        contrasts_to_correct = {'con_0001.nii,1'; 'con_0002.nii,1'; 'con_0003.nii,1'; 'con_0004.nii,1'; 'con_0005.nii,1'; 'con_0006.nii,1';  ...
            'con_0007.nii,1'; 'con_0008.nii,1'; 'con_0009.nii,1'; 'con_0010.nii,1'; 'con_0011.nii,1'; 'con_0012.nii,1'; }
        
        
    for iidx = 1:length(contrasts_to_correct)
        
        clear matlabbatch
        
        
        matlabbatch{1}.spm.util.imcalc.input = {
            fullfile(rootdir, subjects{i}, 'stats/graded_ambiguity/Conventional', contrasts_to_correct{iidx})
            fullfile(rootdir, subjects{i}, 'stats/control/Conventional', contrasts_to_correct{iidx})
            };
        matlabbatch{1}.spm.util.imcalc.output = ['corr' contrasts_to_correct{iidx}];
        matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(rootdir, subjects{i}, 'stats/graded_ambiguity/Conventional')};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1-i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman(run_type,matlabbatch);
    end
    end
    
end

  if define_ANOVA
        
        group_folder = fullfile(rootdir, 'ANOVA');
        
        for i = 1:length(subjects)
        matlabbatch{1}.spm.stats.factorial_design.dir = {group_folder};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Congruency';
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Disparity';
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 6;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0001.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1 2];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0002.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [1 3];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0003.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [1 4];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0004.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).levels = [1 5];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(5).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0005.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).levels = [1 6]
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(6).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0006.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).levels = [2 1];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(7).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0007.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).levels = [2 2];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(8).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0008.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(9).levels = [2 3];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(9).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0009.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(10).levels = [2 4];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(10).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0010.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(11).levels = [2 5];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(11).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0011.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(12).levels = [2 6];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(12).scans(i,1) = {fullfile(rootdir, subjects{i}, 'stats', 'graded_ambiguity', 'Conventional', 'con_0012.nii,1')};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        end
        spm_jobman(run_type,matlabbatch);
    end
    
    
    if estimate_ANOVA
        
        group_folder = fullfile(rootdir, 'ANOVA');
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(group_folder, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        spm_jobman(run_type,matlabbatch);
        
    end

