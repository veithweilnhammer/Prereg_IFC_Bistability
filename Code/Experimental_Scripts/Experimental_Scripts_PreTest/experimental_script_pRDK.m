% Copyright (C) 2019 Veith Weilnhammer, CCM

clear all
close all

%% Pre-Test Experiment

root_dir = '/home/veithweilnhammer/Public/Git/Prereg_IFC_Bistability/Code/Experimental_Scripts/Experimental_Scripts_PreTest/'; % root directory for settings and results
BlueValue = '116' % Put in Blue Value from Heterochromatic flicker photometry
Results.BlueValue = BlueValue;
ObserverName= 'observer_1_veith' % Name of the observer
which_run = 3;12

% R0: dummy run for trying different settings; 
% R1: ambiguous run; 
% R2: fully disambiguated run;
% R3: Percetual Certainty;

if which_run < 3
SettingsName= 'pRDK_pretest'; % points to Settings file
else
SettingsName= 'pRDK_certainty'; % points to Settings file
end

%%
%%
%% Run 0 (Test + Tryout):

if which_run == 0
    session= ['run_' num2str(which_run)]
    run_idx = which_run;
    
    if exist ([root_dir  'Results/Results_' ObserverName '_Conv.mat'], 'file')
        load([root_dir  'Results/Results_' ObserverName '_Conv.mat'])
    end
    
    transition_probability = 0.25; % probability of transition at each overlap. NaN for ambiguous stimulation.
    
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
    
    figure(), set(gcf, 'Color', 'w')
    bar(round(Results.TrialEndTime{1}-Results.TrialStartTime{1})./((mean(mean(frequency(:,[1]))).*Results.n_overlaps*Results.rot_per_trial)+1)), hold on
    ylim([0 40]), xlim([0 2]), ylabel('perceptual phase duration')
    plot([0 2], [35 35], 'r')
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

    
    figure(), set(gcf, 'Color', 'w')
    bar(mean(correct(:,2))), hold on
    ylim([0 1.1]), xlim([0 2]), ylabel('fraction correct')
    plot([0 2], [0.75 0.75], 'r')
    
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
    
    transition_probability = mean(mean(frequency(:,[1]))); % average transition frequency computed in run R1;
    
    clear Results
    [Results] = presentation_pRKD_certainty(root_dir, SettingsName, ObserverName, BlueValue, session, transition_probability);
    
    %% analysis
    [frequency_rating correct_rating rating] = get_conventional_data_rating(Results);
    
    figure(), set(gcf, 'Color', 'w')
    subplot(3,1,1)
    bar(rating(:,1)), ylabel('Uncertainty'), xlim([0 7]), xlabel('Congruent + Level of disambiguation'), ylim([0 4])
    subplot(3,1,2)
    bar(rating(:,2)), ylabel('Uncertainty'), xlim([0 7]), xlabel('Incongruent + Level of disambiguation'), ylim([0 4])
    subplot(3,1,3)
    bar(correct_rating), ylim([0 1.1]), xlim([0 7]), ylabel('fraction correct'), xlabel('Level of disambiguation')
  
    
end



