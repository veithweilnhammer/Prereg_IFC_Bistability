% Copyright (C) 2019 Veith Weilnhammer, CCM
clear all
close all

%% fMRI Experiment

root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Code/Experimental_Scripts/Experimental_Scripts_fMRI/'; % root directory for settings and results
BlueValue = '116' % Put in Blue Value from Heterochromatic flicker photometry
Results.BlueValue = BlueValue;
ObserverName= 'observer_1_veith' % Name of the observer
SettingsName= 'pRDK'; % points to Settings file
which_run = 4; %R0: dummy run for trying different settings; R1: ambiguous runs; R2-R4: runs with parametric disambiguation


%%
%%
%% Run 0 (Test + Tryout):

if which_run == 0
    session= ['run_' num2str(which_run)]
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = 0.5; % probability of transition at each overlap. NaN for ambiguous stimulation.
    
    clear Results
    [Results] = presentation_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability); % run presentation
    
    
    %% analysis
    [test_frequency test_correct] = get_conventional_data(Results); % compute frequency of perceptual events and "correct" perceptual decisions.
end
%%
%%
%% Run 1:
if which_run == 1
    session= ['run_' num2str(which_run)];
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = NaN;
    
    clear Results
    [Results] = presentation_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability);
    
    %% analysis
    [frequency(:,run_idx) correct(:,run_idx)] = get_conventional_data(Results);
    save([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'frequency', 'correct')
    
end

%%
%%
%% Run 2:

if which_run == 2
    session= ['run_' num2str(which_run)];
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = mean(mean(frequency(:,[1]))); % average transition frequency computed in run R1
    
    clear Results
    [Results] = presentation_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability);
    
    %% analysis
    [frequency(:,run_idx) correct(:,run_idx)] = get_conventional_data(Results);
    save([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'frequency', 'correct')
    
    
end
%%
%%
%% Run 3:

if which_run == 3
    session= ['run_' num2str(which_run)];
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = mean(mean(frequency(:,[1])));  % average transition frequency computed in run R1
    
    clear Results
    [Results] = presentation_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability);
    
    %% analysis
    [frequency(:,run_idx) correct(:,run_idx)] = get_conventional_data(Results);
    save([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'frequency', 'correct')
    
    
end
%%
%%
%% Run 4:

if which_run == 4
    session= ['run_' num2str(which_run)];
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_' num2str(session) '.mat'], 'file')
        'Results file already exists. Please change ObserverName or Session'
        return;
    end
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = mean(mean(frequency(:,[1])));  % average transition frequency computed in run R1
    
    clear Results
    [Results] = presentation_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability);
    
    %% analysis
    [frequency(:,run_idx) correct(:,run_idx)] = get_conventional_data(Results);
    save([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'frequency', 'correct')
    
    
end
