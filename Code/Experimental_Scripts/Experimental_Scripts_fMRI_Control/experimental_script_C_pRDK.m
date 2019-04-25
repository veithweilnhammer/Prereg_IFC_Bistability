% Copyright (C) 2019 Veith Weilnhammer, CCM
clear all
close all

%% Control fMRI Experiment

root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Code/Experimental_Scripts/Experimental_Scripts_fMRI_Control/'; % root directory for settings and results
root_dir_main = '/home/veithweilnhammer/Public/Git/Prereg_IFC_BIstability/Code/Experimental_Scripts/Experimental_Scripts_fMRI/' % directory for main experiment (used to load balance betwen perceptual outcomes

BlueValue = '116' % Put in Blue Value from Heterochromatic flicker photometry
Results.BlueValue = BlueValue;
ObserverName= 'observer_1_veith' % Name of the observer
SettingsName= 'C_pRDK'; % points to Settings file
which_run = 1; %R0: dummy run for trying different settings; CR1 control run


%%
%%
%% Run 0 (Test + Tryout):

if which_run == 0
    session= ['run_' num2str(which_run)]
    run_idx = which_run;
    
    %% load data from main experiment
    load(fullfile(root_dir_main, 'Results', ['Results_' ObserverName '_Conv.mat']))
    transition_probability = mean(mean(frequency(:,[1])))
    balance = [mean(correct(:,2:end),2)]' % temporal balance between perceptual interpretations. 0.5: complete balance; 1: complete imbalance;
    
    [Results] = presentation_C_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, balance, transition_probability); % run presentation
    [control_frequency control_correct] = get_conventional_control_data(Results); 
end

%%
%%
%% Control Run 1:
if which_run == 1
    session= ['run_' num2str(which_run)]
    run_idx = which_run;
    
    %% load data from main experiment
    load(fullfile(root_dir_main, 'Results', ['Results_' ObserverName '_Conv.mat']))
    transition_probability = mean(mean(frequency(:,[1])))
    balance = [mean(correct(:,2:end),2)]' % temporal balance between perceptual interpretations. 0.5: complete balance; 1: complete imbalance;
    [Results] = presentation_C_pRKD(root_dir, SettingsName, ObserverName, BlueValue, session, balance, transition_probability); % run presentation
  
    %% analysis
    [control_frequency(:,run_idx) control_correct(:,run_idx)] = get_conventional_control_data(Results);
    save([root_dir  'Results/Results_' ObserverName '_Control_Conv.mat'], 'control_frequency', 'control_correct');
    
end